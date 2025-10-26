from flask import Flask, request, g
from flask_sqlalchemy import SQLAlchemy
from os import path
from flask_login import LoginManager
import os
from websiteNPMINE.config import Config, DevelopmentConfig, ProductionConfig
from flask_bcrypt import Bcrypt
from flask_migrate import Migrate
from flask_bootstrap import Bootstrap4
from flask_mail import Mail
from flask_wtf.csrf import CSRFProtect
import logging
from logging.handlers import RotatingFileHandler
import traceback

db = SQLAlchemy()
migrate = Migrate()
bcrypt = Bcrypt()
login_manager = LoginManager()
login_manager.login_view = 'main.home'
login_manager.login_message_category = 'info'
bootstrap = Bootstrap4()
mail = Mail()
csrf = CSRFProtect()

def create_app(config_class=None):
    app = Flask(__name__)

    if not os.path.exists('logs'):
        os.makedirs('logs')

    logging.basicConfig(
        filename='logs/error.log',
        level=logging.ERROR,
        format='%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
    )

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.ERROR)
    console_handler.setFormatter(logging.Formatter(
        '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
    ))

    @app.teardown_request
    def log_exceptions(error=None):
        """Runs after each request, logs unhandled exceptions."""
        if error is not None:
            app.logger.error(
                "Unhandled Exception: %s\nURL: %s\nMethod: %s\nIP: %s\nTraceback:\n%s",
                error,
                request.url,
                request.method,
                request.remote_addr,
                traceback.format_exc()
            )

    app.logger.addHandler(console_handler)

    if config_class is None:
        env = os.environ.get('FLASK_ENV')
        if env == 'development':
            config_class = DevelopmentConfig
        elif env == 'production':
            config_class = ProductionConfig
        else:
            config_class = DevelopmentConfig 

    app.config.from_object(config_class)
    
    bcrypt.init_app(app)
    db.init_app(app)
    bootstrap.init_app(app)
    mail.init_app(app)
    csrf.init_app(app)
    login_manager.init_app(app)

    from websiteNPMINE.users.routes import users
    from websiteNPMINE.main.routes import main
    from websiteNPMINE.compounds.routes import compounds
    app.register_blueprint(users)
    app.register_blueprint(main)
    app.register_blueprint(compounds)

    migrate.init_app(app, db)
    
    return app


