"""
Django settings for myPhyloDB project.
Generated by 'django-admin startproject' using Django 1.9.1.
For more information on this file, see
https://docs.djangoproject.com/en/1.9/topics/settings/
For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.9/ref/settings/
"""
import os
from django.core.management.utils import get_random_secret_key

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# SECURITY WARNING: keep the secret key used in production secret!
# generate secret key on bootup, will break current session (ie log out all users)
# SECRET_KEY = "ALPHABETSOUPFORTESTINGONLY"   # TODO switch back to random (line below) before commits
SECRET_KEY = get_random_secret_key()   # the actually secure version, will reset logins on boot
ALLOWED_HOSTS = []


# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = False


# Application definition
INSTALLED_APPS = [
    'database',
    'django.contrib.sites',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.admin',
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    #'allauth.socialaccount.providers.google',
]


MIDDLEWARE_CLASSES = [

    # django middleware
    'django.middleware.security.SecurityMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',

    # custom middleware to logout users old sessions when switching machines
    'database.middleware.UserRestrictMiddleware',
]


SITE_ID = 1


ROOT_URLCONF = 'myPhyloDB.urls'


AUTHENTICATION_BACKENDS = (
    # Needed to login by username in Django admin, regardless of `allauth`
    'django.contrib.auth.backends.ModelBackend',
    # `allauth` specific authentication methods, such as login by e-mail
    'allauth.account.auth_backends.AuthenticationBackend',
)


TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            # allauth templates: you could copy this directory into your
            # project and tweak it according to your needs
            # example project specific templates
            os.path.join(BASE_DIR, 'templates'),
            os.path.join(BASE_DIR, 'templates', 'allauth'),
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                # needed for admin templates
                'django.contrib.auth.context_processors.auth',

                # these *may* not be needed
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.contrib.messages.context_processors.messages',

                # allauth needs this from django
                'django.template.context_processors.request',

                # myPhyloDB specific context processors
                'database.views.usrFiles',
            ],
            'libraries': {
                'static': 'django.templatetags.static',
                'i18n': 'django.templatetags.i18n',
                'l10n': 'django.templatetags.l10n',
                'admin_list': 'django.contrib.admin.templatetags.admin_list',
                'admin_modify': 'django.contrib.admin.templatetags.admin_modify',
                'admin_static': 'django.contrib.admin.templatetags.admin_static',
                'admin_urls': 'django.contrib.admin.templatetags.admin_urls',
                'log': 'django.contrib.admin.templatetags.log',

                # allauth template tags
                'account': 'allauth.account.templatetags.account',
                'socialaccount': 'allauth.socialaccount.templatetags.socialaccount'
            }
        },
    }
]


# Django WSGI
WSGI_APPLICATION = 'myPhyloDB.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.9/ref/settings/#databases
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.Microbe'),
    },
    'picrust': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.PICRUSt'),
    }
}

# Password validation
# https://docs.djangoproject.com/en/1.9/ref/settings/#auth-password-validators
AUTH_PASSWORD_VALIDATORS = []


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.9/howto/static-files/
STATIC_URL = '/myPhyloDB/media/'


try:
    from config.allauth_cfg import *
except ImportError:
    pass

try:
    from config.local_cfg import *
except ImportError:
    pass
