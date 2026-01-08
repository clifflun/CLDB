Carvalho Lab DB code base


Sample metadata file format

| pt_id    | family   | project | is_proband | P2_path                                                  | MD_path                                               |
|----------|----------|---------|------------|----------------------------------------------------------|-------------------------------------------------------|
| sample_1| fam_1  | proj_1 | 1          | path/to/P2_file | path/to/MD_file|
| sample_2| fam_1  | proj_1 | 1          | path/to/P2_file | path/to/MD_file|
| sample_3| fam_1  | proj_1 | 1          | path/to/P2_file | path/to/MD_file|


get SLMSuite to local from [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1734-5)

Note: slmseg.py has been modified to make things work

Usage:
Clone this repository git clone
Download docker image from [here](https://hub.docker.com/r/clifflun/sv_db)
Create a docker container and mount the volume to the path where the repository is cloned
Run the docker image in the above container
