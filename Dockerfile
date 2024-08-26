# Start from a base Ubuntu image
FROM ubuntu:latest

# Update and install necessary dependencies
RUN apt-get update && apt-get install -y \
        build-essential \
        cmake \
        wget \
        g++ \
        git \
        python3 \
        python3-pip \
        python3-venv \
        libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Jellyfish dependencies and libraries
RUN apt-get update && apt-get install -y \
        libjellyfish-2.0-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Jellyfish from source
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz \
    && tar -xzf jellyfish-2.3.0.tar.gz \
    && cd jellyfish-2.3.0 \
    && ./configure \
    && make \
    && make install

# Set up a virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install InSilicoSeq and its dependencies within the virtual environment
RUN pip install InSilicoSeq biopython numpy pysam scipy

# Clone and build Naive_Bayes
RUN git clone https://github.com/EESI/Naive_Bayes.git \
    && cd Naive_Bayes \
    && make \
    && cp NB.run /usr/local/bin/NB.run

# Set the working directory
WORKDIR /app
