from django.urls import path
from blog import views

urlpatterns = [
    path('',views.index,name="index"),
    path('articles/<int:id_article>/', views.view_article),
    path('date', views.date_actuelle),
    path('addition/<int:nombre1>/<int:nombre2>/', views.addition),
    path('predictform',views.predictform,name="predictform"),
    path('ajaxpredictform',views.ajaxpredictform,name="ajaxpredictform"),
    path('searchbyprotein',views.searchbyprotein,name="searchbyprotein"),
    path('ajaxsearchprotein',views.ajaxsearchprotein,name="ajaxsearchprotein"),
    path('searchlistprotein',views.searchlistprotein,name="searchlistprotein"),
    path('ajaxsearchlistprotein',views.ajaxsearchlistprotein,name="ajaxsearchlistprotein"),
    path('searchorganism',views.searchorganism,name="searchorganism"),
    path('ajaxsearchorganism',views.ajaxsearchorganism,name="ajaxsearchorganism"),
    path('predicttextmining',views.predicttextmining,name="predicttextmining"),
    path('loadmodel',views.loadmodel,name="loadmodel"),
    path('displayapi',views.displayapi,name="displayapi"),
    path('displaycolab',views.displaycolab,name="displaycolab"),
    path('displaytutorial',views.displaytutorial,name="displaytutorial"),
    path('ajaxpedicttextmining',views.ajaxpedicttextmining,name="ajaxpedicttextmining"),
]
#path('',views.index,name="index"),