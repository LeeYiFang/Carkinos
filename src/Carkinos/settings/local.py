from .base import *

SECRET_KEY = 'br7jtj4&jwn=09^w*882xfcqm94b_#5h5zmtwwm#wtay+yx_0z'
DEBUG = True

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
    }
}
