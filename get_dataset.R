if (!"BiocManager" %in% installed.packages()[, "Package"]) {
    install.packages("BiocManager")
}
    

packages <- c(
    "dplyr", "tibble", "magrittr"
)

for (pkg in packages) {
    if (!pkg %in% installed.packages()[, "Package"])
        BiocManager::install(pkg)
}

suppressMessages({
    library(XMAS2)
    library(dplyr)
    library(tibble)
    library(magrittr)
})

# args

args <- commandArgs(trailingOnly = TRUE)


prepare_lefse <- function(ps,
                          Class,
                          Class_names,
                          Subclass = NULL,
                          cutoff = 10) {
    
    # ps = amplicon_ps_genus
    # Class = "SampleType"
    # Class_names = c("gut", "tongue")
    # Subclass = NULL
    # cutoff = 10
    
    sam_tab <- phyloseq::sample_data(ps) %>%
        data.frame()
    colnames(sam_tab)[which(colnames(sam_tab) == Class)] <- "CompClass"
    
    if (is.null(Subclass)) {
        sam_tab_final <- sam_tab %>%
            dplyr::select(CompClass) %>%
            tibble::rownames_to_column("TempRowNames") %>%
            dplyr::filter(CompClass %in% Class_names) %>%
            dplyr::select(all_of(c("TempRowNames", "CompClass"))) %>%
            tibble::column_to_rownames("TempRowNames")
    } else {
        sam_tab_final <- sam_tab %>%
            dplyr::select(all_of(c("CompClass", Subclass))) %>%
            tibble::rownames_to_column("TempRowNames") %>%
            dplyr::filter(CompClass %in% Class_names) %>%
            dplyr::select(all_of(c("TempRowNames", "CompClass", Subclass))) %>%
            tibble::column_to_rownames("TempRowNames")
    }
    
    colnames(sam_tab_final)[which(colnames(sam_tab_final) == "CompClass")] <- Class
    
    phyloseq::sample_data(ps) <- phyloseq::sample_data(sam_tab_final)
    otu_tab <- phyloseq::otu_table(ps) %>%
        data.frame()
    otu_tab_final <- otu_tab[rowSums(otu_tab) > cutoff, colSums(otu_tab) > cutoff, F]
    phyloseq::otu_table(ps) <- phyloseq::otu_table(as.matrix(otu_tab_final), taxa_are_rows = TRUE)
    
    lefse_data_nosub <- sam_tab_final %>% 
        tibble::rownames_to_column("Sample") %>%
        dplyr::inner_join(otu_tab_final %>% 
                              t() %>% data.frame() %>%
                              tibble::rownames_to_column("Sample"),
                          by = "Sample") %>%
        dplyr::select(-Sample) %>%
        dplyr::select(all_of(Class), all_of(Subclass), everything()) %>%
        t() %>% data.frame()
    
    res <- list(ps=ps,
                lefse=lefse_data,
                lefse_nosub=lefse_data_nosub)
    
    return(res)
}

if (args[1] == "amplicon_ps_genus_lefse") {
    data("amplicon_ps")
    amplicon_ps_genus <- summarize_taxa(amplicon_ps, taxa_level = "Genus")
    amplicon_ps_genus_lefse <- prepare_lefse(
                            ps = amplicon_ps_genus,
                            Class = "SampleType",
                            Class_names = c("gut", "tongue"),
                            cutoff = 10)

    write.table(amplicon_ps_genus_lefse$lefse, 
        "amplicon_ps_genus_lefse.tsv", quote = F, sep = "\t", col.names = F)
    write.table(amplicon_ps_genus_lefse$lefse_nosub, 
        "amplicon_ps_genus_lefse_nosub.tsv", quote = F, sep = "\t", col.names = F)    
} else if (args[1] == "Zeybel_ps_species_lefse") {
    data("Zeybel_Gut")
    Zeybel_ps_species <- summarize_taxa(Zeybel_Gut, taxa_level = "Species")
    Zeybel_ps_species_lefse <- prepare_lefse(
                            ps = Zeybel_ps_species,
                            Class = "LiverFatClass",
                            Class_names = c("Mild", "Moderate"),
                            cutoff = 1e-4)

    write.table(Zeybel_ps_species_lefse$lefse,
        "Zeybel_ps_species_lefse.tsv", quote = F, sep = "\t", col.names = F)
    write.table(Zeybel_ps_species_lefse$lefse_nosub, 
        "Zeybel_ps_species_lefse_nosub.tsv", quote = F, sep = "\t", col.names = F)
}
