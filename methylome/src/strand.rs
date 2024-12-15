#[derive(Hash, PartialEq, Eq, Clone, Copy)]
pub enum Strand {
    Positive,
    Negative,
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let txt = match self {
            Strand::Positive => "+",
            Strand::Negative => "-",
        };
        write!(f, "{}", txt)
    }
}
