import reporter
import dash_auth
import keys

app = reporter.app

dash_auth.BasicAuth(
    app,
    keys.USERNAME_PASSWORD
)

server = app.server