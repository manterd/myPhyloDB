import os

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = '99pbnftb@8#l5$$43)c=c$^y!j4-kz-0t6!zgcp2#k30#q7&g8'

# Application definition
INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.sites',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'registration',
    'database',
    'stats',
    'config',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    #'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    #'django.middleware.security.SecurityMiddleware',
)

ROOT_URLCONF = 'phyloDB.urls'

TEMPLATE_CONTEXT_PROCESSORS = (
    'django.core.context_processors.request',
    'django.contrib.messages.context_processors.messages',
    'django.contrib.auth.context_processors.auth',
)

TEMPLATE_DIRS = (
    os.path.join(BASE_DIR,  'templates'),
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
LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'America/Denver'

USE_I18N = False
USE_L10N = True
USE_TZ = True

# Static files (CSS, JavaScript, Images)
STATIC_URL = '/media/'

STATICFILES_DIRS = (
    os.path.join(BASE_DIR, 'media'),
)
