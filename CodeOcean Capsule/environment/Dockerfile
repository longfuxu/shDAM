# hash:sha256:ceb7d07b76ed34a52de62c233635bff7f1cff682a0d36cbcc714404abcf9d290
FROM registry.codeocean.com/codeocean/miniconda3:4.12.0-python3.9-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN conda install -y \
        jupyter==1.1.1 \
        matplotlib==3.3.4 \
        more-itertools==8.7.0 \
        pandas==1.2.3 \
        scipy==1.12.0 \
        tifffile==2021.3.5 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    jupyterlab==3.0.12 \
    lumicks-pylake==0.8.1 \
    nptdms==1.1.0 \
    opencv-python==4.5.1.48 \
    openpyxl==3.1.5 \
    pwlf==2.2.1 \
    scipy==1.12.0
    
RUN apt-get update && apt-get install -y libgl1-mesa-glx
