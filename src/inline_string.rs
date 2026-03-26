use std::fmt;

const LEN: usize = 24usize;

#[derive(Copy, Clone, PartialEq, Eq, Hash, Ord, PartialOrd)]
#[repr(align(8))]
pub struct InlineString {
    data: [u8; LEN],
}

impl InlineString {
    pub fn new(s: &[u8]) -> Self {
        assert!(
            s.len() <= LEN,
            "The length of the string \"{}\" must be less than or equal to {LEN}",
            std::str::from_utf8(s).unwrap()
        );

        let mut data = [0u8; LEN];
        s.iter().enumerate().for_each(|(i, &b)| data[i] = b);

        Self { data }
    }

    pub fn bytes<'a>(&'a self) -> impl Iterator<Item = u8> + 'a {
        self.data[..self.len()].iter().cloned()
    }

    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.data[..self.len()]).unwrap()
    }

    pub fn len(&self) -> usize {
        let mut len = 0;
        while len < LEN && self.data[len] != 0 {
            len += 1;
        }
        len
    }
}

impl fmt::Debug for InlineString {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\"{}\"",
            std::str::from_utf8(&self.data[..self.len()]).unwrap()
        )
    }
}

impl fmt::Display for InlineString {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            std::str::from_utf8(&self.data[..self.len()]).unwrap()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_empty() {
        let s = InlineString::new(b"");
        assert_eq!(s.len(), 0);
        assert_eq!(s.as_str(), "");
    }

    #[test]
    fn test_new_short() {
        let s = InlineString::new(b"hello");
        assert_eq!(s.len(), 5);
        assert_eq!(s.as_str(), "hello");
    }

    #[test]
    fn test_new_max_len() {
        let data = [b'a'; 24];
        let s = InlineString::new(&data);
        assert_eq!(s.len(), 24);
    }

    #[test]
    #[should_panic]
    fn test_new_too_long() {
        let data = [b'a'; 25];
        InlineString::new(&data);
    }

    #[test]
    fn test_bytes_iterator() {
        let s = InlineString::new(b"ACGT");
        let bytes: Vec<u8> = s.bytes().collect();
        assert_eq!(bytes, b"ACGT");
    }

    #[test]
    fn test_equality() {
        let a = InlineString::new(b"test");
        let b = InlineString::new(b"test");
        let c = InlineString::new(b"other");
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn test_clone_copy() {
        let a = InlineString::new(b"test");
        let b = a;
        assert_eq!(a, b);
    }

    #[test]
    fn test_debug_display() {
        let s = InlineString::new(b"test");
        assert_eq!(format!("{}", s), "test");
        assert_eq!(format!("{:?}", s), "\"test\"");
    }

    #[test]
    fn test_hash() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        set.insert(InlineString::new(b"a"));
        set.insert(InlineString::new(b"b"));
        set.insert(InlineString::new(b"a"));
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn test_ord() {
        let a = InlineString::new(b"aaa");
        let b = InlineString::new(b"bbb");
        assert!(a < b);
    }
}
