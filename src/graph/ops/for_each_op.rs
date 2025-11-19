use crate::graph::*;

pub struct ForEachOp<F: Fn(&mut Read) + Send + Sync> {
    func: F,
}

impl<F: Fn(&mut Read) + Send + Sync> ForEachOp<F> {
    const NAME: &'static str = "ForEachOp";

    /// Apply an arbitrary function on each read.
    pub fn new(func: F) -> Self {
        Self { func }
    }
}

impl<F: Fn(&mut Read) + Send + Sync, T: Trace> GraphNode<T> for ForEachOp<F> {
    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        for read in &mut reads {
            (self.func)(read);
        }
        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}

pub struct DbgOp;

impl DbgOp {
    /// Print each read to standard error.
    pub fn new() -> ForEachOp<impl Fn(&mut Read) + Send + Sync> {
        ForEachOp::new(|read| eprintln!("{read}"))
    }
}

pub struct RemoveInternalOp;

impl RemoveInternalOp {
    /// Remove mappings with labels that start with `_` ("internal" mappings).
    pub fn new() -> ForEachOp<impl Fn(&mut Read) + Send + Sync> {
        ForEachOp::new(|read| read.remove_internal())
    }
}
