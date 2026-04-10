use anyhow::{Result, bail};
use std::{fmt, str::FromStr};

/// Represents a DNA base modification type.
///
/// This enum defines the types of modifications that can occur on DNA bases,
/// including their associated codes for parsing and visualization.
///
/// # Variants
/// - `SixMA`: N6-methyladenine (6mA), represented by the pileup code `a`.
/// - `FiveMC`: 5-methylcytosine (5mC), represented by the pileup code `m`.
/// - `FourMC`: 4-methylcytosine (4mC), represented by the pileup code `21839`.
///
/// # Examples
/// ```
/// use epimetheus_methylome::ModType;
///
/// let mod_type = ModType::SixMA;
/// assert_eq!(mod_type.to_pileup_code(), "a");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash, Copy, PartialOrd, Ord)]
pub enum ModType {
    SixMA,
    FiveMC,
    FourMC,
    FiveHmC,
}

impl ModType {
    /// Returns the pileup code corresponding to the modification type.
    ///
    /// Pileup codes are compact representations of modification types used
    /// in sequencing data (or maybe just modkit):
    /// - `SixMA` (6mA): `"a"`
    /// - `FiveMC` (5mC): `"m"`
    /// - `FourMC` (4mC): `"21839"`
    ///
    /// # Examples
    /// ```
    /// use epimetheus_methylome::ModType;
    ///
    /// let mod_type = ModType::FiveMC;
    /// assert_eq!(mod_type.to_pileup_code(), "m");
    /// ```
    pub fn to_pileup_code(&self) -> &'static str {
        match self {
            ModType::SixMA => "a",
            ModType::FiveMC => "m",
            ModType::FourMC => "21839",
            ModType::FiveHmC => "h",
        }
    }

    /// Passes the sam header tag to ModType
    ///
    /// Sam tags have the following and more for the methylation types:
    /// MM:Z:{*}
    /// - `A+a` (6mA): `"a"`
    /// - `C+m` (5mC): `"m"`
    /// - `C+21839` (4mC): `"21839"`
    ///
    /// # Examples
    /// ```
    /// use epimetheus_methylome::ModType;
    ///
    /// assert_eq!(ModType::from_sam_code('A', "a"), Some(ModType::SixMA));
    /// ```
    pub fn from_sam_code(base: char, modification: &str) -> Option<Self> {
        // WARN I should add the other bases in the future
        match (base, modification) {
            ('A', "a") => Some(ModType::SixMA),
            ('C', "m") => Some(ModType::FiveMC),
            ('C', "21839") => Some(ModType::FourMC),
            _ => None,
        }
    }
}

impl fmt::Display for ModType {
    /// Formats the modification type for display purposes.
    ///
    /// Each modification type is represented in the format:
    /// `<Modification Name> (<Pileup Code>)`.
    ///
    /// For example:
    /// - `6mA (a)` for `SixMA`
    /// - `5mC (m)` for `FiveMC`
    /// - `4mC (21839)` for `FourMC`
    ///
    /// # Examples
    /// ```
    /// use epimetheus_methylome::ModType;
    ///
    /// let mod_type = ModType::FourMC;
    /// assert_eq!(format!("{}", mod_type), "4mC (21839)");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ModType::SixMA => write!(f, "6mA (a)"),
            ModType::FiveMC => write!(f, "5mC (m)"),
            ModType::FourMC => write!(f, "4mC (21839)"),
            ModType::FiveHmC => write!(f, "5hmC (h)"),
        }
    }
}

/// Parses a modification type from a string.
///
/// The input string must match one of the following:
/// - `"a"` for `SixMA` (6mA)
/// - `"m"` for `FiveMC` (5mC)
/// - `"21839"` for `FourMC` (4mC)
///
/// # Arguments
/// - `mod_type`: A string slice representing the modification type.
///
/// # Returns
/// - `Ok(ModType)` if the string matches a supported modification type.
/// - `Err` if the string does not match any supported modification type.
///
/// # Examples
/// ```
/// use epimetheus_methylome::ModType;
///
/// let mod_type = "a".parse::<ModType>().unwrap();
/// assert_eq!(mod_type, ModType::SixMA);
///
/// let invalid = "unsupported".parse::<ModType>();
/// assert!(invalid.is_err());
/// ```
impl FromStr for ModType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "a" => Ok(ModType::SixMA),
            "m" => Ok(ModType::FiveMC),
            "21839" => Ok(ModType::FourMC),
            "h" => Ok(ModType::FiveHmC),
            _ => bail!("Unsupported mod type: {}", s),
        }
    }
}
