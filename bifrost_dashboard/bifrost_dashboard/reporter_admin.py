import os
import dash_auth
import bifrost_dashboard.reporter
import yaml
config = yaml.safe_load(open(os.environ["BIFROST_DASH_CONFIG"]))

bifrost_dashboard.reporter.ADMIN = True

app = bifrost_dashboard.reporter.app

dash_auth.BasicAuth(
    app,
    config['USERNAME_PASSWORD']
)

server = app.server
