use std::fmt::Display;

use anyhow::bail;

/// Represents an IUPAC nucleotide base.
///
/// IUPAC nucleotide codes are used to represent ambiguous positions in DNA or RNA sequences.
/// This enum includes both standard nucleotides (A, T, G, C) and ambiguous codes (e.g., R, Y, N).
///
/// # Variants
/// - `A`: Adenine
/// - `T`: Thymine
/// - `G`: Guanine
/// - `C`: Cytosine
/// - `R`: Purine (A or G)
/// - `Y`: Pyrimidine (C or T)
/// - `S`: Strong interaction (C or G)
/// - `W`: Weak interaction (A or T)
/// - `K`: Keto (G or T)
/// - `M`: Amino (A or C)
/// - `B`: Not A (C, G, or T)
/// - `D`: Not C (A, G, or T)
/// - `H`: Not G (A, C, or T)
/// - `V`: Not T (A, C, or G)
/// - `N`: Any base
///
/// # References
/// Based on IUPAC nucleotide code conventions.
/// For more details, see: https://en.wikipedia.org/wiki/Nucleic_acid_notation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IupacBase {
    A,
    T,
    G,
    C,
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
    N,
}

impl Display for IupacBase {
    /// Formats the `IupacBase` as a string for display purposes.
    ///
    /// This implementation converts the enum variant to its corresponding IUPAC character.
    ///
    /// # Example
    /// ```
    /// use methylome::IupacBase;
    ///
    /// let base = IupacBase::A;
    /// assert_eq!(format!("{}", base), "A");
    /// ```
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

impl IupacBase {
    /// Parses a single character into an `IupacBase` enum variant.
    ///
    /// # Arguments
    /// - `base`: A character representing an IUPAC nucleotide code.
    ///
    /// # Returns
    /// - `Ok(IupacBase)` if the character is a valid IUPAC nucleotide code.
    /// - `Err` if the character is invalid.
    ///
    /// # Examples
    /// ```
    /// use methylome::IupacBase;
    ///
    /// let base = IupacBase::parse_char('A').unwrap();
    /// assert_eq!(base, IupacBase::A);
    ///
    /// let invalid = IupacBase::parse_char('Z');
    /// assert!(invalid.is_err());
    /// ```
    pub fn parse_char(base: char) -> anyhow::Result<Self> {
        let iupac_base = match base {
            'A' => Ok(Self::A),
            'T' => Ok(Self::T),
            'G' => Ok(Self::G),
            'C' => Ok(Self::C),
            'R' => Ok(Self::R),
            'Y' => Ok(Self::Y),
            'S' => Ok(Self::S),
            'W' => Ok(Self::W),
            'K' => Ok(Self::K),
            'M' => Ok(Self::M),
            'B' => Ok(Self::B),
            'D' => Ok(Self::D),
            'H' => Ok(Self::H),
            'V' => Ok(Self::V),
            'N' => Ok(Self::N),
            _ => bail!("Not a defined Iupac base: {base}"),
        };
        iupac_base
    }

    /// Converts the `IupacBase` to its string representation.
    ///
    /// This method maps the `IupacBase` variant to its corresponding single-character string.
    ///
    /// # Examples
    /// ```
    /// use methylome::IupacBase;
    ///
    /// let base = IupacBase::A;
    /// assert_eq!(base.to_string(), "A");
    /// ```
    pub fn to_string(&self) -> String {
        match self {
            IupacBase::A => "A".to_string(),
            IupacBase::T => "T".to_string(),
            IupacBase::G => "G".to_string(),
            IupacBase::C => "C".to_string(),
            IupacBase::R => "R".to_string(),
            IupacBase::Y => "Y".to_string(),
            IupacBase::S => "S".to_string(),
            IupacBase::W => "W".to_string(),
            IupacBase::K => "K".to_string(),
            IupacBase::M => "M".to_string(),
            IupacBase::B => "B".to_string(),
            IupacBase::D => "D".to_string(),
            IupacBase::H => "H".to_string(),
            IupacBase::V => "V".to_string(),
            IupacBase::N => "N".to_string(),
        }
    }

    /// Returns the complement of the given IUPAC base.
    ///
    /// Complements are defined as follows:
    /// - A ↔ T
    /// - G ↔ C
    /// - R ↔ Y
    /// - S ↔ S
    /// - W ↔ W
    /// - K ↔ M
    /// - B ↔ V
    /// - D ↔ H
    /// - N ↔ N
    ///
    /// # Arguments
    /// - `base`: A reference to an `IupacBase`.
    ///
    /// # Returns
    /// The complement of the input base as an `IupacBase`.
    ///
    /// # Examples
    /// ```
    /// use methylome::IupacBase;
    ///
    /// let complement = IupacBase::to_complement_base(&IupacBase::A);
    /// assert_eq!(complement, IupacBase::T);
    /// ```
    pub fn to_complement_base(base: &IupacBase) -> Self {
        match base {
            IupacBase::A => IupacBase::T,
            IupacBase::T => IupacBase::A,
            IupacBase::G => IupacBase::C,
            IupacBase::C => IupacBase::G,
            IupacBase::R => IupacBase::Y,
            IupacBase::Y => IupacBase::R,
            IupacBase::S => IupacBase::S,
            IupacBase::W => IupacBase::W,
            IupacBase::K => IupacBase::M,
            IupacBase::M => IupacBase::K,
            IupacBase::B => IupacBase::V,
            IupacBase::D => IupacBase::H,
            IupacBase::H => IupacBase::D,
            IupacBase::V => IupacBase::B,
            IupacBase::N => IupacBase::N,
        }
    }

    /// Converts the `IupacBase` into its corresponding regular expression representation.
    ///
    /// This is useful for pattern matching ambiguous nucleotide sequences.
    ///
    /// # Examples
    /// ```
    /// use methylome::IupacBase;
    ///
    /// let regex = IupacBase::R.to_regex();
    /// assert_eq!(regex, "[AG]");
    /// ```
    pub fn to_regex(&self) -> &str {
        match self {
            IupacBase::A => "A",
            IupacBase::T => "T",
            IupacBase::G => "G",
            IupacBase::C => "C",
            IupacBase::R => "[AG]",
            IupacBase::Y => "[CT]",
            IupacBase::S => "[CG]",
            IupacBase::W => "[AT]",
            IupacBase::K => "[GT]",
            IupacBase::M => "[AC]",
            IupacBase::B => "[CGT]",
            IupacBase::D => "[AGT]",
            IupacBase::H => "[ACT]",
            IupacBase::V => "[ACG]",
            IupacBase::N => ".",
        }
    }
}
