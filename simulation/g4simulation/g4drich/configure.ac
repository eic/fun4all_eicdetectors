AC_INIT(g4drich,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)

LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
dnl if test $ac_cv_prog_gxx = yes; then
dnl    CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-unused-variable -Wno-unused-but-set-variable"
dnl fi
case $CXX in
  clang++)
   CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-macro-redefined"
  ;;
  *g++)
   CXXFLAGS="$CXXFLAGS -Wall -Werror"
  ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
