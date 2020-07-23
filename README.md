# RNA-seq-2020

This repo contains all the code used to analyze and create figures for the 2020 RNA-seq analysis on PC and RACM.

## Format

There are several folders to organize specific code:

* DE-analysis - contains all code used for edgeR/deseq2 analysis
* heatmaps - contains all code used for heatmap creation
* venn - contains all code used for venn diagram creation
* bar-graphs - contains all code used for bar graph creation

## Contributing

There are two ways to add files to the repo:

1. Publish straight to the website

In Github, go to the appropriate folder and upload the file you want to share. When you add it, include a short description of the file before uploading it.

2. Publish from terminal

I prefer uploading the files straight from terminal, since working in terminal makes it easier to upload and download the files. To do this....

* First install [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
* Add your github account info with the following commands `git config --global user.name "John Doe"` and `git config --global user.email johndoe@example.com`, swapping John Doe and johndoe@example.com with your github username and associated email
* **Clone** or download the repository using the following command `git clone https://github.com/vedantham-lab/RNA-seq-2020.git`. This command will create a new folder in whatever directory you are currently in in terminal, so make sure you are `cd`'d into the directory you want the folder to be inside of.
* Move or create code inside one of the appropriate folders of the repo
* **Add** your new files. You can add all the files at once with `git add -A` or add each individual file with `git add file1.txt file2.txt`. This tells git that these are files it should monitor for changes.
* **Commit** your files. This will tell git that these are files you want to upload. `git commit -m "write a short description here of your files!`"
* **Push** your files. This will upload your files straight to github so everyone else can view it. `git push origin master`

And to download new files from the repo

* **Pull** the files to your computer. `git pull origin master`

## Guidelines

Git can get confused if multiple people upload the same file with different contents. So when uploading your code or figures, add your initials to the end of the file.

Also, don't upload any significant data to github (no readcount tables, fastq, bam files) 

