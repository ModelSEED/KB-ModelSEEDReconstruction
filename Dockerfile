FROM kbase/sdkpython:3.8.0
MAINTAINER chenry@anl.gov
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install --upgrade pip

# update security libraries in the base image
RUN pip install cffi --upgrade \
    && pip install pyopenssl --upgrade \
    && pip install ndg-httpsclient --upgrade \
    && pip install pyasn1 --upgrade \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade \
	&& pip install --upgrade pyopenssl ndg-httpsclient && \
    pip install --upgrade pyasn1 requests 'requests[security]' && \
    pip install coverage networkx cython && \
    pip install --upgrade pip setuptools wheel cffi

RUN apt-get update
RUN apt-get install -y gcc
RUN rm -rf /miniconda/lib/python3.6/site-packages/numpy
RUN rm -rf /miniconda/lib/python3.6/site-packages/ruamel*
RUN pip install --upgrade pip
RUN pip install "numpy<1.24"
RUN pip install cobra
RUN pip install networkx
RUN pip install chemw==0.3.2

RUN echo '13' >/dev/null && pip install --use-deprecated=legacy-resolver git+https://github.com/cshenry/ModelSEEDpy.git
RUN echo '7' >/dev/null && pip install git+https://github.com/cshenry/cobrakbase.git@5c0bfefb569a2540df878fdf995889590412f232

RUN echo '16' >/dev/null && mkdir deps && cd deps && \
	git clone --branch main https://github.com/cshenry/KBBaseModules.git
RUN mkdir test

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]