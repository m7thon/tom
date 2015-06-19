from distutils.core import setup

setup(
    name = "tom",
    version = "0.4.0",
    author = "Michael Thon",
    author_email = "m7.thon@gmail.com",
    description = ("toolkit for observable operator modeling"),
    license = "MIT",
    url = "https://gitlab.com/m7.thon/tom",
    packages=['tom'],
    package_data={'tom': ['_tomlib.so'],},
)
