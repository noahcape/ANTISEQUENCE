use cfg_if::cfg_if;

pub trait SeedSearcher {
    fn search(&self, text: &[u8], candidate_fn: impl FnMut(SeedMatch));
}

pub struct SmallSearcher<const K: usize> {
    pattern_luts: [Aligned<64>; K],
    pattern_idxs: Vec<(u32, u32)>,
}

impl SmallSearcher<const K: usize> {
    pub fn new(patterns: impl Iterator<Item = (usize, &[u8])>) -> Result<Self, ()> {
        cfg_if! {
            if #[cfg(not(feature = "avx2"))] {
                return Err(());
            }
        }

        let mut pattern_luts = [Aligned::<64>::default(); K];
        let mut pattern_idxs = Vec::new();
        let mut idx = 0;

        for (&pattern_idx, pattern) in patterns.iter() {
            for (pattern_i, kmer) in pattern.windows(K).enumerate() {
                if idx >= 8 {
                    return Err(());
                }

                for (i, &c) in kmer.iter().enumerate() {
                    let c = if i % 2 == 0 { (c as usize) & 0b0000_1111 } else { ((c as usize) >> 1) & 0b0000_1111 };

                    for j in 0..4 {
                        pattern_luts[i].0[c + j * 16] |= (1 << idx) as u8;
                    }
                }

                pattern_idxs.push((pattern_idx, pattern_i));
                idx += 1;
            }
        }

        Ok(Self { pattern_luts, pattern_idxs })
    }
}

impl<const K: usize> SeedSearcher for SmallSearcher<{ K }> {
    fn search(&self, text: &[u8], mut candidate_fn: impl FnMut(SeedMatch)) {
        cfg_if! {
            if #[cfg(feature = "avx2")] {
                const L: usize = 32;

                #[inline(always)]
                fn inner<const REST: bool>(text: &[u8], candidate_fn: &mut impl FnMut(SeedMatch), i: usize, len: usize) {
                    let lo_4_mask = _mm256_set1_epi8(0b0000_1111i8);
                    let mut set = _mm256_setzero_si256();

                    let load_mask = if REST {
                        load_mask_avx2(len)
                    } else {
                        _mm256_set1_epi8(-1i8)
                    };

                    for j in 0..K {
                        let mut chars = if REST {
                            load_rest_avx2(text.as_ptr().add(i + j), len)
                        } else {
                            _mm256_loadu_si256(text.as_ptr().add(i + j) as _)
                        };

                        if j % 2 != 0 {
                            chars = _mm256_srli_epi16(chars, 1);
                        }

                        let pattern_lut = _mm256_load_si256(self.pattern_luts[j].0.as_ptr() as _);
                        let curr_set = _mm256_shuffle_epi8(pattern_lut, _mm256_and_si256(chars, lo_4_mask));
                        set = _mm256_and_si256(set, curr_set);
                    }

                    set = _mm256_and_si256(set, load_mask);
                    let nonzero = _mm256_cmpgt_epi8(set, _mm256_setzero_si256());
                    let mut nonzero_mask = _mm256_movemask_epi8(nonzero) as u32;

                    if nonzero_mask > 0 {
                        let mut a = Aligned::<{ L }>::default();
                        _mm256_store_si256(a.0.as_mut_ptr() as _, set);

                        while nonzero_mask > 0 {
                            let idx = nonzero_mask.trailing_zeros() as usize;
                            let s = *a.0.as_ptr().add(idx) as usize;
                            let hash_idx = s.trailing_zeros() as usize;
                            let (pattern_idx, pattern_i) = *self.pattern_idxs.as_ptr().add(hash_idx);
                            candidate_fn(SeedMatch { pattern_idx: pattern_idx as usize, pattern_i: pattern_i as usize, text_i: i + idx });

                            nonzero_mask &= nonzero_mask - 1;
                        }
                    }
                }

                let mut i = 0;

                while i + (K - 1) + L <= text.len() {
                    inner::<false>(text, &mut candidate_fn, i, L);
                    i += L;
                }

                if i + (K - 1) < text.len() {
                    inner::<true>(text, &mut candidate_fn, i, text.len() - i - (K - 1));
                }
            } else {
                unreachable!()
            }
        }
    }
}

pub struct GeneralSearcher {
    k: usize,
    table: HashToPatternIdx,
}

impl GeneralSearcher {
    const B: usize = 8;

    pub fn new(patterns: impl Iterator<Item = (usize, &[u8])>, k: usize) -> Self {
        let mut hashes = Vec::new();

        for (&pattern_idx, p) in patterns.iter() {
            Self::get_hashes(p, k, |curr_hashes, len, pattern_i| {
                hashes.extend(curr_hashes[..len].into_iter().enumerate().map(|(i, h)| (h, pattern_idx, pattern_i + i)));
            });
        }

        Self {
            k,
            table: HashToPatternIdx::new(hashes, (hashes.len() * 4).next_power_of_two()),
        }
    }

    #[inline(always)]
    unsafe fn get_hashes(s: &[u8], k: usize, f: impl FnMut([u32; Self::B], usize, usize)) {
        if s.len() < k {
            return;
        }

        let mut hash = 0;
        let mut i = 0;
        let ptr = s.as_ptr();

        while i < k {
            hash = hash.rotate_left(1) ^ wyhash_byte(*ptr.add(i));
            i += 1;
        }

        while i + Self::B <= s.len() {
            let mut hashes = [0u32; Self::B];
            for j in 0..Self::B {
                hashes[j] = hash;
                hash = hash.rotate_left(1) ^ wyhash_byte(*ptr.add(i + j)) ^ wyhash_byte(*ptr.add(i + j - k)).rotate_left(k);
            }

            f(hashes, Self::B, i);

            i += B;
        }

        let len = s.len() - i;

        if len > 0 {
            let mut hashes = [0u32; Self::B];
            for j in 0..len {
                hashes[j] = hash;
                hash = hash.rotate_left(1) ^ wyhash_byte(*ptr.add(i + j)) ^ wyhash_byte(*ptr.add(i + j - k)).rotate_left(k);
            }

            f(hashes, len, i);
        }
    }
}

impl SeedSearcher for GeneralSearcher {
    fn search(&self, text: &[u8], candidate_fn: impl FnMut(SeedMatch)) {
        Self::get_hashes(text, self.k, |curr_hashes, len, text_i| {
            self.table.batch_lookup::<{ Self::B }>(curr_hashes, len, text_i, &mut candidate_fn);
        });
    }
}

struct HashToPatternIdx {
    hash_to_idxs: Vec<u32>,
    pattern_idxs: Vec<(u32, u32)>,
}

impl HashToPatternIdx {
    pub fn new(mut pattern_idxs: Vec<(u32, usize, usize)>) -> Self {
        let num_hashes = (pattern_idxs.len() * 8).next_power_of_two();
        pattern_idxs.sort_unstable();
        let mut hash_to_idxs = vec![0u32; num_hashes];

        for (hash, _, _) in &pattern_idxs {
            hash_to_idxs[(hash as usize) & (hash_to_idxs.len() - 1)] += 1;
        }

        let mut sum = 0;
        hash_to_idxs.iter_mut().for_each(|i| {
            let temp = *i;

            if temp == 0 {
                *i = std::u32::MAX;
            } else {
                *i = sum;
                sum += temp;
            }
        });

        let pattern_idxs = pattern_idxs.into_iter().map(|(_, pattern_idx, pattern_i)| (pattern_idx, pattern_i)).collect();

        Self {
            hash_to_idxs,
            pattern_idxs,
        }
    }

    #[inline(always)]
    pub unsafe fn batch_lookup<const N: usize>(mut hashes: [usize; N], len: usize, text_i: usize, candidate_fn: &mut impl FnMut(SeedMatch)) {
        let mut starts = [0u32; N];
        let mut found = 0;

        for i in 0..N {
            hashes[i] = if i < len { hashes[i] & (self.hash_to_idxs.len() - 1) } else { 0 };
            starts[i] = *self.hash_to_idxs.as_ptr().add(hashes[i]);
            found += if starts[i] == std::u32::MAX { 0 } else { 1 };
        }

        if found == 0 {
            return;
        }

        for (i, (start, hash)) in starts[..len].into_iter().zip(hashes[..len].into_iter()).enumerate() {
            if start == std::u32::MAX {
                continue;
            }

            let start = start as usize;
            let end = self.hash_to_idxs.get(hash + 1).map(|i| i as usize).unwrap_or(self.pattern_idxs.len());

            for j in start..end {
                let (pattern_idx, pattern_i) = *self.pattern_idxs.as_ptr().add(j);
                candidate_fn(SeedMatch { text_i: text_i + i, pattern_idx: pattern_idx as usize, pattern_i: pattern_i as usize });
            }
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct SeedMatch {
    pub pattern_idx: usize,
    pub pattern_i: usize,
    pub text_i: usize,
}

#[derive(Default)]
#[align(64)]
struct Aligned<const L: usize>([u8; L]);

cfg_if! {
    if #[cfg(feature = "avx2")] {
        fn load_mask_avx2(len: usize) -> __m256i {
            static MASKS: [u8; 64] = {
                let mut m = [0u8; 64];
                let mut i = 0;

                while i < 32 {
                    m[i] = 0xFFu8;
                    i += 1;
                }

                m
            };
            _mm256_loadu_si256(MASKS.as_ptr().add(32 - len) as _)
        }

        fn read_rest_avx2(ptr: *const u8, len: usize) -> __m256i {
            let addr = ptr as usize;
            let start_page = addr >> 12;
            let end_page = (addr + 31) >> 12;

            if start_page == end_page {
                _mm256_loadu_si256(ptr as _)
            } else {
                let mut a = Aligned::<32>([0u8; 32]);
                std::ptr::copy_nonoverlapping(a.0.as_mut_ptr(), ptr, len);
                _mm256_load_si256(a.0.as_ptr() as _)
            }
        }
    }
}

#[inline(always)]
fn wyhash_byte(b: u8) -> u32 {
    let a = 0xa076_1d64_78bd_642fu64;
    let b = (b as u64) ^ 0xe703_7ed1_a0b4_28dbu64;
    let c = a.wrapping_mul(b);
    (c ^ (c >> 32)) as u32
}
