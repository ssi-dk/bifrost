FROM ssidk/bifrost-base:2.0.5

ARG version="2.0.5"
ARG last_updated="08/06/2020"
ARG name="run_launcher"
ARG full_name="bifrost-${name}"

LABEL \
    name=${name} \
    description="Docker environment for ${full_name}" \
    version=${version} \
    resource_version=${last_updated} \
    maintainer="kimn@ssi.dk;"

#- Tools to install:start---------------------------------------------------------------------------
# None
#- Tools to install:end ----------------------------------------------------------------------------

#- Additional resources (files/DBs): start ---------------------------------------------------------
# None
#- Additional resources (files/DBs): end -----------------------------------------------------------

#- Source code:start -------------------------------------------------------------------------------
RUN cd /bifrost && \
    git clone --branch ${version} https://github.com/ssi-dk/${full_name}.git ${name};
#- Source code:end ---------------------------------------------------------------------------------

#- Set up entry point:start ------------------------------------------------------------------------
ENV PATH /bifrost/${name}/:$PATH
ENTRYPOINT ["launcher.py"]
CMD ["launcher.py", "--help"]
#- Set up entry point:end --------------------------------------------------------------------------