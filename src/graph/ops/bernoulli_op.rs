use rand::distributions::Bernoulli;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::graph::*;

pub struct BernoulliOp {
    attr: Attr,
    bernoulli: Bernoulli,
    seed: u64,
}

impl BernoulliOp {
    const NAME: &'static str = "BernoulliOp";

    /// Set the attribute `attr` to a sampled boolean from a Bernoulli distribution
    /// with probability `prob` of true.
    ///
    /// This is fully deterministic for a chosen seed and ordering of reads, even with
    /// multiple threads.
    pub fn new(attr: Attr, prob: f64, seed: u32) -> Self {
        Self {
            attr,
            bernoulli: Bernoulli::new(prob)
                .unwrap_or_else(|e| panic!("Error creating bernoulli distribution: {e}")),
            seed: seed as u64,
        }
    }
}

impl<T: Trace> GraphNode<T> for BernoulliOp {
    fn run_inner(&self, mut read: Read) -> Result<(Option<Read>, bool)> {
        // use the index of the read in the seed for determinism when multithreading
        let seed = (self.seed << 32).wrapping_add(read.first_idx() as u64);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        let rand_bool = self.bernoulli.sample(&mut rng);

        // panic to make borrow checker happy
        *read
            .data_mut(self.attr.str_type, self.attr.label, self.attr.attr)
            .unwrap_or_else(|e| panic!("Error in {}: {e}", Self::NAME)) = Data::Bool(rand_bool);

        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn finish(&self) -> Result<bool> {
        Ok(true)
    }
}
