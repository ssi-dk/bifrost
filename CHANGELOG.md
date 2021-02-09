# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.2.0] - 2020-10-07
### Added
- CHANGELOG.md into repo

### Changed
- Templated files to be common for SSI maintained pipelines. All files use .env and <COMPONENT_NAME>/config.yaml as the source of all information. The config.yaml should be considered the primary source of all information regarding the component and it's settings. The .env file needs to contain the <COMPONENT_NAME> and install specific settings (currently just mongo_db connection). 
  - docker-compose.dev.yaml
  - docker-compose.yaml
  - .env
    - This is being used for both Dockerfile and passing the values into the Docker image, not sure if that has any issues.
  - setup.py
    - This can't use libraries to extract config values, so right now it's hardcoded on what to look for. This can cause some potential issues.
- The following files are also impacted by the changes
  - Dockerfile

### Removed
- Docker-compose files no longer point to a custom env file