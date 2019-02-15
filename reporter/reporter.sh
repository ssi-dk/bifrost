. activate bifrost
uwsgi -s /tmp/reporter.sock --processes 2 --manage-script-name --mount /=reporter:server --http :8050