from forms import UserRegForm
from models import UserProfile


def user_created(sender, user, request, **kwargs):
    form = UserRegForm(request.POST)
    data = UserProfile(user=user)
    data.firstName = form.data['firstName'] or 'No data'
    data.lastName = form.data['lastName'] or 'No data'
    data.affiliation = form.data['affiliation'] or 'No data'
    data.city = form.data['city'] or 'No data'
    data.state = form.data['state'] or 'No data'
    data.country = form.data['country'] or 'No data'
    data.zip = form.data['zip'] or 'No data'
    data.phone = form.data['phone'] or 'No data'
    data.reference = form.data['reference'] or 'No data'
    data.purpose = form.data['purpose'] or 'No data'
    data.save()

from registration.signals import user_registered
user_registered.connect(user_created)
