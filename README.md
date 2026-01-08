# Description

CLDB is a Python-based toolkit designed to support flexible management and querying of copy-number variant (CNV) / structural variant (SV) metadata for genomic analyses, particularly in the context of rare disease cohorts. It provides utilities for ingesting sample metadata, filtering and querying CNV/SV call sets, and integrating with visualization and analytic workflows.

This codebase includes modules for querying, helper utilities, web pages, and miscellaneous scripts to build reproducible analysis pipelines.

# Usage

1. Clone the Repository

```bash
git clone https://github.com/clifflun/CLDB.git
cd CLDB
```

2. Download the Docker Image

Download the prebuilt Docker image from the provided source (see Docker Image section below).

```bash
docker pull clifflun/sv_db
```

3. Run the Docker Container

Start a Docker container and mount the local CLDB repository into the container. Expose port 8501 to access the web application:

```
docker run -it \
  -p 8501:8501 \
  -v /path/to/CLDB:/CLDB \
 clifflun/sv_db 
```

4. Access the Application

Once the container is running, the application will start automatically.
Open a web browser and navigate to:
```
http://localhost:8501
```

