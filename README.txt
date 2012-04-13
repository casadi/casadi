On x86_64 platforms there may be problem linking problems when you compile Sundials on your own.
To solve this, configure Sundials with the following options:

./configure --with-cflags="-O3 -fPIC"


