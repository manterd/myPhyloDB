import os

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = '99pbnftb@8#l5$$43)c=c$^y!j4-kz-0t6!zgcp2#k30#q7&g8'

# Application definition
INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'registration',
    'south',
    'phyloDB',
    'config',
    'database',
    'database.anova',
    'database.diffabund',
    'database.export',
    'database.pcoa',
    'database.spls',
    'stats',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'phyloDB.urls'

TEMPLATE_CONTEXT_PROCESSORS = (
    'django.core.context_processors.request',
    'django.contrib.messages.context_processors.messages',
    'django.contrib.auth.context_processors.auth',
)

TEMPLATE_DIRS = (
    os.path.join(BASE_DIR, 'templates'),
)

SESSION_EXPIRE_AT_BROWSER_CLOSE = True

WSGI_APPLICATION = 'phyloDB.wsgi.application'

# Database
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'dbMicrobe')
    },
}

# Internationalization
USE_L10N = True

USE_I18N = True
LANGUAGE_CODE = 'en-us'

USE_TZ = True
TIME_ZONE = 'America/Denver'

# Static files (CSS, JavaScript, Images)
STATIC_ROOT = BASE_DIR
STATIC_URL = '/media/'

STATICFILES_DIRS = (
    os.path.join(BASE_DIR, 'media'),
)
