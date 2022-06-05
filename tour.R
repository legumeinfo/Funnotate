# ----------------------------------------------------------
# library(cicerone)
# ----------------------------------------------------------

tour <- Cicerone$new()$
  step("tour-intro", "Phylogram Tour",
    "This quick tour will acquaint you with the phylogeny tree viewer and other resources available in this section. Use the Next button or → (right arrow key) to advance the tour. Use the Prev button or ← (left arrow key) to step back.")$
  step("phylogramFamilyInfo", "Gene Family Name",
    "The gene family name and description are shown in the header. Family descriptions are derived from homology-based functional analysis of the hidden Markov model representing the known sequences in the family, and include InterPro and Gene Ontology identifiers.")$
  step("tour-phylotree", "Phylogram",
    "The phylogram view displays a phylogenetic tree with branch spans proportional to the amount of character change, including both legumes and selected non-legumes.")$
  step(".leaf.tnt_tree_node", "Terminal Nodes",
    "The tree nodes on the right are terminal nodes. Click on a (colored) node to see more information about legume genes including links to organism, gene, genomic context, and various resources and viewers. For example, the genomic context viewer shows flanking genes and allows exploration of syntenic regions from all included legume genomes.", position = "right-center", is_id = FALSE)$
  step(".inner.tnt_tree_node", "Interior Nodes",
    "The other tree nodes are interior nodes. Click on a (white) interior node to view genomic context and genomic distribution for the node's sub-tree. In addition, you can expand and collapse subtrees and you can focus the phylogram to the selected subtree. If you expand, collapse, or focus on a subtree, the changes are instantly reflected in the taxon view and in the MSA view (more about those in a moment).", position = "right-center", is_id = FALSE)$
  step(".root.tnt_tree_node", "Root Node",
    "The root node is also an interior node, although it is colored black for reference. (Note: this root is only one of several possible root choices, and may not be the oldest common ancestor. It is the result of midpoint rooting of the tree.)", position = "right-center", is_id = FALSE)$
  step("div.shiny-options-group > label:nth-of-type(1)", "Taxa and Legend",
    "The Taxa and Legend view shows the taxonomic distribution of the current set features from the phylogram, and provides a legend for the legumes (colored dots). You can toggle the Taxa and Legend view via this checkbox.", is_id = FALSE)$
  step("tour-taxaAndLegend", "Taxa and Legend view",
    "You can click the species names to toggle them on and off. You can double-click a species name to filter to only that species. Changes are instantly reflected in the Phylogram view and in the MSA view.", position = "top")$
  step("div.shiny-options-group > label:nth-of-type(2)", "MSA Visualization",
    "The Multiple Sequence Alignment for the family is available via this checkbox.", is_id = FALSE)$
  step("tour-msa", "Multiple Sequence Alignment",
    "The MSA view shows the current set of features from the phylogram, plus the consensus sequence. You can select one or more feature names in the MSA, and they will be highlighted in the Phylogram view.", position = "top")$
  step("tour", "Done", "Click Done to end the tour.", position = "right")

# --------------------------------------------------------------

