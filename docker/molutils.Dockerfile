From python:3.9.16-slim-bullseye
RUN pip install --no-cache-dir moltemplate==2.20.19
RUN apt update && apt install wget build-essential gfortran bc gawk openbabel -y && \
    wget https://github.com/m3g/packmol/archive/refs/tags/v20.14.2.tar.gz && \
    tar -xzf v20.14.2.tar.gz && cd packmol-20.14.2 && make && cp packmol /usr/bin/ && \
    cd ../ && rm -r v20.14.2.tar.gz packmol-20.14.2 && \
    apt-get purge -y wget build-essential && apt-get autoremove -y && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
