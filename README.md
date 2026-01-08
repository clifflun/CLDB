Carvalho Lab DB code base


Sample metadata file format

| pt_id    | family   | project | is_proband | P2_path                                                  | MD_path                                               |
|----------|----------|---------|------------|----------------------------------------------------------|-------------------------------------------------------|
| sample_1| fam_1  | proj_1 | 1          | path/to/P2_file | path/to/MD_file|
| sample_2| fam_1  | proj_1 | 1          | path/to/P2_file | path/to/MD_file|
| sample_3| fam_1  | proj_1 | 1          | path/to/P2_file | path/to/MD_file|

## Usage

### 1. Clone the Repository

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
  -v /path/to/CLDB:/app/CLDB \
 clifflun/sv_db 
```

4. Access the Application

Once the container is running, the application will start automatically.
Open a web browser and navigate to:
```
http://localhost:8501
```
