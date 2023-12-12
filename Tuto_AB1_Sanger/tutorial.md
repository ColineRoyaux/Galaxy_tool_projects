---
layout: tutorial_hands_on

title: 'Clean and manage Sanger sequences from raw files to aligned consensus'
zenodo_link: https://zenodo.org/records/7104640/files/AOPEP_and_CHD8_sequences_20220907.zip?download=1
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- contributor1
- contributor2

---


# Introduction

<!-- This is a comment. 

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!
-->
## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) :
>
>    ```
>    https://zenodo.org/records/7104640/files/AOPEP_and_CHD8_sequences_20220907.zip?download=1
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Check that the datatype is '.zip'
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 4. Create primer FASTA file, copy:
>    ```
>    >Forward_CHD8
>    GAGGTGAAAGAATCATAAATTGG
>    >Reverse_CHD8
>    CCCTGTGTACAAATAGCTTTTGT
>    ```
>    - Open the Galaxy Upload Manager ({% icon galaxy-upload %} on the top-right of the tool panel)
>    - Select **Paste/Fetch Data**
>    - Paste into the text field
>    - Change **Type (set all):** from "Auto-detect" to `fasta`
>    - Change the name from "New File" to "Primer file"
>    - Click **Start**
>
{: .hands_on}

# Prepare primer data

<!--It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:
-->

## Separate and format primers files

> <hands-on-title> Create separate files for each primer </hands-on-title>
>
> 1. {% tool [Filter FASTA](toolshed.g2.bx.psu.edu/repos/galaxyp/filter_by_fasta_ids/filter_by_fasta_ids/2.3) %} with the following parameters:
>    - {% icon param-file %} *"FASTA sequences"*: `Primer file`
>    - *"Criteria for filtering on the headers"*: `Regular expression on the headers`
>        - *"Regular expression pattern the header should match"*: `Reverse_CHD8`
>    - Add tags "#Primer" and "#Reverse"
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 2. {% tool [Filter FASTA](toolshed.g2.bx.psu.edu/repos/galaxyp/filter_by_fasta_ids/filter_by_fasta_ids/2.3) %} with the following parameters:
>    - {% icon param-file %} *"FASTA sequences"*: `Primer file`
>    - *"Criteria for filtering on the headers"*: `Regular expression on the headers`
>        - *"Regular expression pattern the header should match"*: `Forward_CHD8`
>    - Add tags "#Primer" and "#Forward"
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 3. Remove eventual gaps from primers {% tool [Degap.seqs](toolshed.g2.bx.psu.edu/repos/iuc/mothur_degap_seqs/mothur_degap_seqs/1.39.5.0) %} with the following parameters:
>    - {% icon param-files %} *"fasta - Dataset"*: `Two Filter FASTA outputs` (outputs of **Filter FASTA** {% icon tool %})
>    {% snippet faqs/galaxy/tools_select_multiple_datasets.md %}
>
{: .hands_on}



> <hands-on-title> Compute Reverse-Complement of the Reverse primer </hands-on-title>
>
> 1. {% tool [Reverse-Complement](toolshed.g2.bx.psu.edu/repos/devteam/fastx_reverse_complement/cshl_fastx_reverse_complement/1.0.2+galaxy0) %} the sequence Reverse primer with the following parameters:
>    - {% icon param-file %} *"Input file in FASTA or FASTQ format"*: `Degap.seqs #Reverse FASTA output` (output of **Degap.seqs** {% icon tool %})
>
{: .hands_on}



# Prepare sequence data

## Unzip data files

> <hands-on-title> Unzip </hands-on-title>
>
> 1. {% tool [Unzip](toolshed.g2.bx.psu.edu/repos/imgteam/unzip/unzip/6.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"input_file"*: `AOPEP_and_CHD8_sequences_20220907.zip?download=1`
>    - *"Extract single file"*: `All files`
>
> > ### {% icon question %} Question
> > How many files is there in the ZIP archive ?
> >
> > > ### {% icon solution %} Solution
> > > 12 (if you have a different number of files something likely went srong)
> > >
> > {: .solution}
> >
> {: .question}
>
{: .hands_on}

From now on, we'll be working a lot on data collections:
{% snippet faqs/galaxy/tools_select_collection.md %}


## Filter collection to separate forward and reverse sequence files

> <hands-on-title> Filter </hands-on-title>
>
> 1. {% tool [Extract element identifiers](toolshed.g2.bx.psu.edu/repos/iuc/collection_element_identifiers/collection_element_identifiers/0.0.2) %} with the following parameters:
>    - {% icon param-collection %} *"Dataset collection*: `output collection` (output of **Unzip** {% icon tool %})
>
> 2. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output` (output of **Extract element identifiers** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^[A-Za-z0-9_-]+F$`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^[A-Za-z0-9_-]+AOPEP[A-Za-z0-9_-]+$`
>    - Tag output with "#Reverse"
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 3. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output` (output of **Extract element identifiers** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^[A-Za-z0-9_-]+R$`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^[A-Za-z0-9_-]+AOPEP[A-Za-z0-9_-]+$`
>    - Tag output with "#Forward"
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 4. {% tool [Filter collection](__FILTER_FROM_FILE__) %} with the following parameters:
>    - {% icon param-collection %} *"Input Collection*: `output collection` (output of **Unzip** {% icon tool %})
>    - *"How should the elements to remove be determined?"*: `Remove if identifiers are ABSENT from file`
>        - {% icon param-file %} *"Filter out identifiers absent from"*: `#Forward files list` (output of **Regex Find And Replace** {% icon tool %})
>
> 5. {% tool [Filter collection](__FILTER_FROM_FILE__) %} with the following parameters:
>    - {% icon param-collection %} *"Input Collection*: `output collection` (output of **Unzip** {% icon tool %})
>    - *"How should the elements to remove be determined?"*: `Remove if identifiers are ABSENT from file`
>        - {% icon param-file %} *"Filter out identifiers absent from"*: `#Reverse files list` (output of **Regex Find And Replace** {% icon tool %})
>
>    > <comment-title> What's happening in this section? </comment-title>
>    >
>    > First step: Extracting the list of file names in the data collection
>    > Second step: Removing file names with a "F" and "AOPEP" in their name -> creating a list of Reverse sequence files of the marker CH8
>    > Third step: Removing file names with a "R" and "AOPEP" in their name -> creating a list of Forward sequence files of the marker CH8
>    > Fourth and fifth step: Select files in the collection -> creating two distinct collections with Forward sequence files on one hand and Reverse sequence file on the other hand
>    > 
>    {: .comment}
>
{: .hands_on}

## Convert AB1 sequence files to FASTQ and trim low-quality ends

> <hands-on-title> AB1 to FASTQ files and trim low quality ends </hands-on-title>
> 
> Do these steps twice !! We have Froward and Reverse sequence data collections, do these steps starting with each "(filtered)" data collections
>
> 1. {% tool [ab1 to FASTQ converter](toolshed.g2.bx.psu.edu/repos/ecology/ab1_fastq_converter/ab1_fastq_converter/1.20.0) %} with the following parameters:
>    - {% icon param-collection %} *"Input ab1 file"*: `(filtered) output collection` (output of **Filter collection** {% icon tool %})
>    - *"Do you want trim ends according to quality scores ?"*: `No, use full sequences.`
>
> 2. {% tool [seqtk_trimfq](toolshed.g2.bx.psu.edu/repos/iuc/seqtk/seqtk_trimfq/1.3.1) %} with the following parameters:
>    - {% icon param-collection %} *"Input FASTA/Q file"*: `output collection` (output of **ab1 to FASTQ converter** {% icon tool %})
>    - *"Mode for trimming FASTQ File"*: `Quality`
>        - *"Maximally trim down to INT bp"*: `0`
>
{: .hands_on}

## Compute reverse complement sequence for Reverse sequences only 

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTQ Groomer](toolshed.g2.bx.psu.edu/repos/devteam/fastq_groomer/fastq_groomer/1.1.5) %} with the following parameters:
>    - {% icon param-collection %} *"File to groom"*: `output collection` (output of **seqtk_trimfq** {% icon tool %})
>    - *"Advanced Options"*: `Show Advanced Options`
>        - *"Summarize input data"*: `Do not Summarize Input (faster)`
>
>    > <comment-title> What is this step? </comment-title>
>    >
>    > It is just a necessary step to get the right input format for the following step **Reverse-Complement** {% icon tool %}
>    {: .comment}
>
> 2. {% tool [Reverse-Complement](toolshed.g2.bx.psu.edu/repos/devteam/fastx_reverse_complement/cshl_fastx_reverse_complement/1.0.2+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input file in FASTA or FASTQ format"*: `output collection` (output of **FASTQ Groomer** {% icon tool %})
>
>
{: .hands_on}

## Merge corresponding Forward and Reverse sequences single files

> <hands-on-title> Sort collections </hands-on-title>
> 
> Do this step twice !! One has to make sure Forward and Reverse sequences collections are in the same order to get the right forward and the right reverse sequence to be merged together
>
> 1. {% tool [Sort collection](__SORTLIST__) %} with the following parameters:
>    - {% icon param-collection %} *"Input Collection"*: `Collection` (output of **seqtk_trimfq** {% icon tool %} & output of **Reverse-Complement** {% icon tool %})
>    - *"Sort type"*: `alphabetical`
>
{: .hands_on}

> <hands-on-title> Merge Forward and Reverse sequence files </hands-on-title>
>
> 1. {% tool [seqtk_mergepe](toolshed.g2.bx.psu.edu/repos/iuc/seqtk/seqtk_mergepe/1.3.1) %} with the following parameters:
>    - {% icon param-file %} *"Input FASTA/Q file #1"*: `output` (output of **Sort collection** {% icon tool %})
>    - {% icon param-file %} *"Input FASTA/Q file #2"*: `output` (output of **Sort collection** {% icon tool %})
>
> Check there is two sequences in each three files of the newly-created collection
>
{: .hands_on}

## Convert FASTQ files to FASTA

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTQ Groomer](toolshed.g2.bx.psu.edu/repos/devteam/fastq_groomer/fastq_groomer/1.1.5) %} with the following parameters:
>    - {% icon param-collection %} *"File to groom"*: `default` (output of **seqtk_mergepe** {% icon tool %})
>    - *"Advanced Options"*: `Show Advanced Options`
>        - *"Summarize input data"*: `Do not Summarize Input (faster)`
>
>    > <comment-title> What is this step? </comment-title>
>    >
>    > It is just a necessary step to get the right input format for the following step **FASTQ to FASTA** {% icon tool %}
>    {: .comment}
>
> 1. {% tool [FASTQ to FASTA](toolshed.g2.bx.psu.edu/repos/devteam/fastq_to_tabular/fastq_to_tabular/1.1.5) %} with the following parameters:
>    - {% icon param-collection %} *"FASTQ file to convert"*: `output collection` (output of **FASTQ Groomer** {% icon tool %})
>    - *"Discard sequences with unknown (N) bases"*: `no`
>    - *"Rename sequence names in output file (reduces file size)"*: `no`
>    - *"Compress output FASTA"*: `No`
>
>    > <comment-title> information </comment-title>
>    >
>    > If this step doesn't work, one can try tools **FASTQ to tabular** {% icon tool %} and **tabular to FASTA** {% icon tool %} instead
>    {: .comment}
>
{: .hands_on}

## Align sequences and retrieve consensus for each sequence

> <hands-on-title> Align and consensus </hands-on-title>
>
> 1. {% tool [Align sequences](toolshed.g2.bx.psu.edu/repos/iuc/qiime_align_seqs/qiime_align_seqs/1.9.1.0) %} with the following parameters:
>    - {% icon param-collection %} *"Input fasta file"*: `output` (output of **Tabular-to-FASTA** {% icon tool %})
>    - *"Method for aligning sequences"*: `clustalw`
>    - *"Minimum percent sequence identity to closest blast hit to include sequence in alignment"*: `0.1`
>
> 2. {% tool [Consensus sequence from aligned FASTA](toolshed.g2.bx.psu.edu/repos/ecology/aligned_to_consensus/aligned_to_consensus/1.0.0) %} with the following parameters:
>    - {% icon param-collection %} *"Input fasta file with at least two sequences"*: `aligned_sequences` (output of **Align sequences** {% icon tool %})
>
> 3. {% tool [Merge.files](toolshed.g2.bx.psu.edu/repos/iuc/mothur_merge_files/mothur_merge_files/1.39.5.0) %} with the following parameters:
>    - *"Merge"*: `fasta files`
>        - {% icon param-collection %} *"inputs - fasta"*: `output` (output of **Consensus sequence from aligned FASTA** {% icon tool %})
>
{: .hands_on}

# Primers AND sequences 

## Merge and align consensus sequences file and primer files

> <hands-on-title> Merge and format consensus sequences + primers file </hands-on-title>
>
> 1. {% tool [Merge.files](toolshed.g2.bx.psu.edu/repos/iuc/mothur_merge_files/mothur_merge_files/1.39.5.0) %} with the following parameters:
>    - *"Merge"*: `fasta files`
>        - {% icon param-files %} *"inputs - fasta"*: `consensus sequences` (output of **Merge.files** {% icon tool %}), `Reverse primer` (output of **Reverse-Complement** {% icon tool %}), `Forward primer` (output of **Degap.seqs** {% icon tool %})
>
> 2. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output` (output of **Merge.files** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `([A-Z-])>`
>            - *"Replacement"*: `\1\n>`
>
>    > <comment-title> What's going on in this second step? </comment-title>
>    >
>    > Sometimes, **Merge.files** {% icon tool %} doesn't keep line feed between the sequences, this step permits to correct it and get a FASTA file that is formatted properly 
>    {: .comment}
>
{: .hands_on}

> <hands-on-title> Align sequences and primers </hands-on-title>
>
> 1. {% tool [Align sequences](toolshed.g2.bx.psu.edu/repos/iuc/qiime_align_seqs/qiime_align_seqs/1.9.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input fasta file"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Method for aligning sequences"*: `mafft`
>    - *"Minimum percent sequence identity to closest blast hit to include sequence in alignment"*: `0.1`
>
{: .hands_on}

## Check your sequencing hasn't been contaminated and your sequence belongs to the right group by computing a BLAST on the NCBI database

> <hands-on-title> NVBI Blast </hands-on-title>
>
> 1. {% tool [NCBI BLAST+ blastn](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastn_wrapper/2.10.1+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Nucleotide query sequence(s)"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"Subject database/sequences"*: `Locally installed BLAST database`
>        - *"Nucleotide BLAST database"*: ``
>    - *"Output format"*: `Tabular (select which columns)`
>        - *"Standard columns"*: ``
>        - *"Extended columns"*: ``
>        - *"Other identifier columns"*: ``
>    - *"Advanced Options"*: `Show Advanced Options`
>        - *"Maximum hits to consider/show"*: `10`
>        - *"Restrict search of database to a given set of ID's"*: `No restriction, search the entire database`
>
> > ### {% icon question %} Question
> > The sequences we cleaned belong to what species?
> >
> > > ### {% icon solution %} Solution
> > > *Homo sapiens*
> > >
> > {: .solution}
> >
> {: .question}
>
{: .hands_on}

# Conclusion

We successfully cleaned AB1 sequence files !