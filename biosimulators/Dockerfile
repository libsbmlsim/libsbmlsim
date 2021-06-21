# Base OS
FROM python:3.9-slim-buster

ARG VERSION=1.4.0
ARG LIBSBML_VERSION=5.19.0

# metadata
LABEL \
    org.opencontainers.image.title="LibSBMLsim" \
    org.opencontainers.image.version="${VERSION}" \
    org.opencontainers.image.description="Library for simulating an SBML model which contains Ordinary Differential Equations" \
    org.opencontainers.image.url="http://fun.bio.keio.ac.jp/software/libsbmlsim/" \
    org.opencontainers.image.documentation="http://fun.bio.keio.ac.jp/software/libsbmlsim/" \
    org.opencontainers.image.source="https://github.com/libsbmlsim/libsbmlsim" \
    org.opencontainers.image.authors="LibSBMLsim Team <info@fun.bio.keio.ac.jp>" \
    org.opencontainers.image.vendor="LibSBMLsim Team" \
    org.opencontainers.image.licenses="LGPL-2.1-or-later" \
    \
    base_image="python:3.9-slim-buster" \
    version="${VERSION}" \
    software="LibSBMLsim" \
    software.version="${VERSION}" \
    about.summary="Library for simulating an SBML model which contains Ordinary Differential Equations" \
    about.home="http://fun.bio.keio.ac.jp/software/libsbmlsim/" \
    about.documentation="http://fun.bio.keio.ac.jp/software/libsbmlsim/" \
    about.license_file="https://github.com/libsbmlsim/libsbmlsim/blob/master/LICENSE.txt" \
    about.license="SPDX:LGPL-2.1-or-later" \
    about.tags="BioSimulators,mathematical model,kinetic model,simulation,systems biology,computational biology,SBML,SED-ML,COMBINE,OMEX" \
    maintainer="LibSBMLsim Team <info@fun.bio.keio.ac.jp>"

# install libSBML
RUN apt-get update -y \
    && apt-get install --no-install-recommends -y \
        wget \
        libxml2 \
    \
    && wget https://master.dl.sourceforge.net/project/sbml/libsbml/${LIBSBML_VERSION}/stable/Linux/64-bit/libSBML-${LIBSBML_VERSION}-Linux-x64.deb \
    && dpkg -i libSBML-${LIBSBML_VERSION}-Linux-x64.deb \
    && rm libSBML-${LIBSBML_VERSION}-Linux-x64.deb \
    && export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH \
    \
    && apt-get remove -y \
        wget \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
ENV LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

# install libSBMLSim
COPY . /root/libsbmlsim/
RUN apt-get update -y \
    && apt-get install --no-install-recommends -y \
        wget \
        build-essential \
        cmake \
        swig \
        libbz2-dev \
    \
    && cd /root/libsbmlsim \
    && mkdir build \
    && cd build \
    && cmake \
        -D LIBSBML_INCLUDE_DIR=/usr/include \
        -D LIBSBML_LIBRARY=/usr/lib64/libsbml.so \
        -D WITH_PYTHON=ON \
        .. \
    && make \
    && make install \
    && cd /root \
    && rm -r libsbmlsim \
    \
    && apt-get remove -y \
        wget \
        build-essential \
        cmake \
        swig \
        libbz2-dev \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Copy code for command-line interface into image and install it
COPY biosimulators/ /root/biosimulators/
RUN pip install /root/biosimulators \
    && rm -rf /root/biosimulators

# supported environment variables
ENV ALGORITHM_SUBSTITUTION_POLICY=SIMILAR_VARIABLES \
    VERBOSE=0 \
    MPLBACKEND=PDF

# Entrypoint
ENTRYPOINT ["biosimulators-libsbmlsim"]
CMD []
