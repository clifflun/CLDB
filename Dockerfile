# Use a Python base image
FROM clifflun/sv_db:base

# RUN apt-get update && apt-get install -y \
#    r-base \
#    r-base-dev \
#    libssl-dev \
#    libxml2-dev \
#    libcurl4-openssl-dev
    
# RUN R -e "install.packages('ggplot2', repos='http://cran.r-project.org')"
# RUN R -e "install.packages('dplyr', repos='http://cran.r-project.org')"


# Set the working directory inside the container
WORKDIR /CLDB

# Copy the requirements.txt file and install dependencies
COPY requirements.txt .

RUN pip install -r requirements.txt

# Copy the Streamlit app files into the container

COPY Hello.py .
COPY .streamlit/ ./.streamlit

# Expose the Streamlit port
EXPOSE 8501

# Command to run Streamlit app
CMD ["streamlit", "run", "Hello.py"]
