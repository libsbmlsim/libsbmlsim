# Base OS
FROM python:2.7-slim-buster

ARG VERSION=1.4.0

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
    base_image="python:2.7-slim-buster" \
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

# Install requirements
RUN wget http://fun.bio.keio.ac.jp/software/libsbmlsim/downloads/libsbmlsim_${VERSION}_amd64.deb -o /tmp/libsbmlsim_${VERSION}_amd64.deb \
    sudo dpkg -i /tmp/libsbmlsim_${VERSION}_amd64.deb \
    rm /tmp/libsbmlsim_${VERSION}_amd64.deb

# Copy code for command-line interface into image and install it
COPY biosimulators/ /root/biosimulators/
RUN pip install /root/biosimulators \
    && rm -rf /root/biosimulators

# supported environment variables
ENV ALGORITHM_SUBSTITUTION_POLICY=SIMILAR_VARIABLES \
    VERBOSE=0

# Entrypoint
ENTRYPOINT ["libsbmlsim"]
CMD []
