use std::collections::HashSet;
use std::convert::From;
use std::fmt::{Display, Formatter};
use std::iter::Iterator;

pub trait IsGap {
    fn is_gap(&self) -> bool;
}
#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub enum Base {
    A,
    C,
    T,
    G,
    GAP,
    F,
    L,
    I,
    M,
    V,
    S,
    P,
    Y,
    STOP,
    H,
    Q,
    N,
    K,
    D,
    E,
    W,
    R,
    UNKNOWN,
}

#[derive(Debug)]
pub struct MultipleSequenceAlignment {
    pub(crate) columns: Vec<Vec<SequenceElement>>,
    pub(crate) num_seqs: usize,
    pub(crate) width: usize,
}

impl From<char> for Base {
    fn from(c: char) -> Base {
        match c {
            'A' => Base::A,
            'C' => Base::C,
            'T' => Base::T,
            'G' => Base::G,
            '-' => Base::GAP,
            'F' => Base::F,
            'L' => Base::L,
            'I' => Base::I,
            'M' => Base::M,
            'V' => Base::V,
            'S' => Base::S,
            'P' => Base::P,
            'Y' => Base::Y,
            '*' => Base::STOP,
            'H' => Base::H,
            'Q' => Base::Q,
            'N' => Base::N,
            'K' => Base::K,
            'D' => Base::D,
            'E' => Base::E,
            'W' => Base::W,
            'R' => Base::R,
            _ => Base::UNKNOWN,
        }
    }
}

impl Into<char> for Base {
    fn into(self) -> char {
        match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::T => 'T',
            Base::G => 'G',
            Base::GAP => '-',
            Base::F => 'F',
            Base::L => 'L',
            Base::I => 'I',
            Base::M => 'M',
            Base::V => 'V',
            Base::S => 'S',
            Base::P => 'P',
            Base::Y => 'Y',
            Base::STOP => '*',
            Base::H => 'H',
            Base::Q => 'Q',
            Base::N => 'N',
            Base::K => 'K',
            Base::D => 'D',
            Base::E => 'E',
            Base::W => 'W',
            Base::R => 'R',
            Base::UNKNOWN => 'X',
        }
    }
}

impl Into<char> for &Base {
    fn into(self) -> char {
        match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::T => 'T',
            Base::G => 'G',
            Base::GAP => '-',
            Base::F => 'F',
            Base::L => 'L',
            Base::I => 'I',
            Base::M => 'M',
            Base::V => 'V',
            Base::S => 'S',
            Base::P => 'P',
            Base::Y => 'Y',
            Base::STOP => '*',
            Base::H => 'H',
            Base::Q => 'Q',
            Base::N => 'N',
            Base::K => 'K',
            Base::D => 'D',
            Base::E => 'E',
            Base::W => 'W',
            Base::R => 'R',
            Base::UNKNOWN => 'X',
        }
    }
}
#[derive(Debug)]
pub struct Sequence {
    pub(crate) elements: Vec<SequenceElement>,
}

impl Iterator for Sequence {
    type Item = SequenceElement;

    fn next(&mut self) -> Option<Self::Item> {
        self.elements.iter().next().cloned()
    }
}

impl FromIterator<SequenceElement> for Sequence {
    fn from_iter<T: IntoIterator<Item = SequenceElement>>(iter: T) -> Self {
        Sequence {
            elements: iter.into_iter().collect(),
        }
    }
}

#[derive(Debug, Hash, Copy, Clone, Eq)]
pub struct SequenceElement {
    sequence_id: Option<usize>,
    base: Base,
    pub(crate) position: Option<usize>,
}

impl PartialEq<Self> for SequenceElement {
    fn eq(&self, other: &Self) -> bool {
        self.base.eq(&other.base)
            && self.position.eq(&other.position)
            && self.sequence_id.eq(&other.sequence_id)
    }
}

impl Display for SequenceElement {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let s_id: String = match self.sequence_id {
            None => String::from("_"),
            Some(c) => c.to_string(),
        };

        let pos: String = match self.position {
            None => String::from("_"),
            Some(c) => c.to_string(),
        };

        write!(f, "[{:?}(s{}-p{})]", self.base, s_id, pos)
    }
}
impl Display for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let c: char = self.into();
        write!(f, "{:?}", c)
    }
}

impl Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let output_str: String = self.elements.iter().map(|e| format!("{}", e)).collect();
        write!(f, "{}", output_str)
    }
}
impl IsGap for SequenceElement {
    fn is_gap(&self) -> bool {
        matches!(self.base, Base::GAP)
    }
}

pub trait ElementFromString {
    fn from_characters(s: Vec<u8>, id: &usize) -> Self;
}

impl ElementFromString for Sequence {
    fn from_characters(s: Vec<u8>, id: &usize) -> Sequence {
        let mut seq: Vec<SequenceElement> = Vec::with_capacity(s.len());
        let mut count: usize = 0;
        for c in s.iter() {
            let base = Base::from(*c as char);
            match base {
                Base::GAP => {
                    seq.push(SequenceElement {
                        sequence_id: Some(*id),
                        base,
                        position: count.checked_sub(1),
                    });
                }
                _ => {
                    seq.push(SequenceElement {
                        sequence_id: Some(*id),
                        base,
                        position: Some(count),
                    });
                    count += 1;
                }
            }
        }

        Sequence { elements: seq }
    }
}

pub type MsaHashSets<'a> = Vec<Vec<Option<HashSet<&'a SequenceElement>>>>;
