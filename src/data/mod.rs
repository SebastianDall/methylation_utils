pub mod contig;
pub mod methylation;

use crate::data::contig::Contig;
use std::collections::HashMap;

pub struct GenomeWorkspace {
    pub contigs: HashMap<String, Contig>,
}

impl GenomeWorkspace {
    pub fn new() -> Self {
        Self {
            contigs: HashMap::new(),
        }
    }

    pub fn add_contig(&mut self, contig: Contig) {
        self.contigs.insert(contig.id.clone(), contig);
    }

    pub fn prune_empty_contigs(&mut self) {
        self.contigs
            .retain(|_id, contig| !contig.methylated_positions.is_empty());
    }

    pub fn get_mut_contig(&mut self, id: &str) -> Option<&mut Contig> {
        self.contigs.get_mut(id)
    }
}
