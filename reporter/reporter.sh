. activate bifrost
uwsgi -s /tmp/reporter.sock --manage-script-name --mount /=reporter:server --http :8050