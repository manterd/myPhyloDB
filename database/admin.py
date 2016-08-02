from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from django.contrib.auth.models import User
from database.models import user_profile


admin.site.unregister(User)


class UserProfileInline(admin.StackedInline):
    model = user_profile


class UserProfileAdmin(UserAdmin):
    inlines = [UserProfileInline, ]


admin.site.register(User, UserProfileAdmin)
