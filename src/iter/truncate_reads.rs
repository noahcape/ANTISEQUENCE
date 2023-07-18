use crate::iter::*;

pub struct TruncateReads<R: Reads> {
    reads: R,
    selector_expr: SelectorExpr,
    labels: Vec<Label>,
    by_length: EndIdx,
}

impl<R: Reads> TruncateReads<R> {
    pub fn new(
        reads: R,
        selector_expr: SelectorExpr,
        labels: Vec<Label>,
        by_length: EndIdx,
    ) -> Self {
        Self {
            reads,
            selector_expr,
            labels,
            by_length,
        }
    }
}

impl<R: Reads> Reads for TruncateReads<R> {
    fn next_chunk(&self) -> Result<Vec<Read>> {
        let mut reads = self.reads.next_chunk()?;

        for read in reads.iter_mut() {
            if !(self
                .selector_expr
                .matches(read)
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: "truncate reads",
                })?)
            {
                continue;
            }

            self.labels
                .iter()
                .try_for_each(|l| read.truncate(l.str_type, l.label, self.by_length))
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: "truncate reads",
                })?;
        }

        Ok(reads)
    }

    fn finish(&mut self) -> Result<()> {
        self.reads.finish()
    }
}
