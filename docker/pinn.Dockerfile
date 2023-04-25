FROM tensorflow/tensorflow:2.8.4-gpu

RUN apt update \
    && apt install pip git -y \
    && pip install \
      git+https://github.com/Teoroo-CMC/PiNN.git@9911d59 \
      git+https://github.com/yqshao/tips.git@3fa9326 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
