use crate::iter::*;

pub struct TruncateToReads<R: Reads> {
    reads: R,
    selector_expr: SelectorExpr,
    labels: Vec<Label>,
    to_length: EndIdx,
}

impl<R: Reads> TruncateToReads<R> {
    pub fn new(
        reads: R,
        selector_expr: SelectorExpr,
        labels: Vec<Label>,
        to_length: EndIdx,
    ) -> Self {
        Self {
            reads,
            selector_expr,
            labels,
            to_length,
        }
    }
}

impl<R: Reads> Reads for TruncateToReads<R> {
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
                .try_for_each(|l| read.truncate_to(l.str_type, l.label, self.to_length))
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
