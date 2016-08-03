from forms import UserRegForm
from models import UserProfile


def user_created(sender, user, request, **kwargs):
    print "user_created:", request.POST
    form = UserRegForm(request.POST)
    data = UserProfile(user=user)
    data.affiliation = form.data['affiliation']
    data.city = form.data['city']
    data.state = form.data['state']
    data.country = form.data['country']
    data.zip = form.data['zip']
    data.phone = form.data['phone']
    data.reference = form.data['reference']
    data.purpose = form.data['purpose']
    data.save()


from registration.signals import user_registered
user_registered.connect(user_created)