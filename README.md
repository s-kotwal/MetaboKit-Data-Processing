# MetaboKit-Data-Processing 240414

MetaboKit (https://github.com/MetaboKit/MetaboKit)

Presented is a code used to automatically filter .tsv files "Report_RTseparated" or "Report_RTseparated_fill" that are auto-generated from processing Mass Spec .mzML files using MetaboKit. The final output is a file that removes duplicate annotations, filter compounds based blanks, abundance cutoffs, etc. 

Requirements: R (4.4.2), MetaboKit (DDA v20240513). 

Refer to 'example_files.zip' for example script and .tsv files.
Must have blank and experimental samples in MetaboKit .tsv file as 'pos-(sample)-it1' to indicate mode, sample name, and number of iterative injection.
