@conf_vars=qw( CXXFLAGS LDFLAGS LIBS );

sub allowed_options {
   my ($allowed_options, $allowed_with)=@_;
   # @$allowed_options{ qw( ACTION ) }=();
   @$allowed_with{ qw( mongoc mongoc-include mongoc-lib bson-include bson-lib ) }=();
}

sub usage {
   print STDERR " --with-mongoc=PATH          installation path of mongoc and bson libraries, if non-standard\n",
                " --with-mongoc-include=PATH  include path of mongoc library\n",
                " --with-mongoc-lib=PATH      path with mongoc library\n",
                " --with-bson-include=PATH    include path of bson library\n",
                " --with-bson-lib=PATH        path with bson library\n",
}


sub proceed {
   my ($options)=@_;
   my $path;
   my $lmc1 = "libmongoc-1.0";
   my $lb1 = "libbson-1.0";
   my $usedopts = 0;

   $CXXFLAGS  = "-Wno-shadow ";

   if (defined($path = $options->{mongoc})){
      $usedopts++;
      if (-d "$path/include/$lmc1" && -d "$path/include/$lb1") {
         $CXXFLAGS  = "-I$path/include/$lmc1";
         $CXXFLAGS .= " -I$path/include/$lb1";
         if ($path ne "/usr") {
            my $libpath=Polymake::Configure::get_libdir($path, "mongoc-1.0");
            $LDFLAGS = "-L$libpath -Wl,-rpath,$libpath";
         }
      } else {
         # make sure we fail the main configure here instead of just disabling
         # the interface
         $options->{polydb} = "";
         die <<"---";
libmongoc or libbson headers not found in $path/include. Make sure 
$path/include/$lmc1
$path/include/$lb1
are valid.
---
      }
   } else {
      if (defined($path = $options->{"mongoc-include"})){
         die "at 2";
         $usedopts++;
         $path .= "/$lmc1" if -d "$path/$lmc1";
         $CXXFLAGS .= " -I$path";
      }
      if (defined($path = $options->{"bson-include"})){
         die "at 3";
         $usedopts++;
         $path .= "/$lb1" if -d "$path/$lb1";
         $CXXFLAGS .= " -I$path";
      }
      if (defined($path = $options->{"mongoc-lib"})){
         die "at 4";
         $usedopts++;
         $LDFLAGS .= " -L$path";
         if ($path ne "/usr") {
            $LDFLAGS .= " -Wl,-rpath,$path";
         }
      }
      if (defined($path = $options->{"bson-lib"})){
         die "at 5";
         $usedopts++;
         $LDFLAGS .= " -L$path";
         if ($path ne "/usr") {
            $LDFLAGS .= " -Wl,-rpath,$path";
         }
      }
   }
   $LIBS = "-lmongoc-1.0 -lbson-1.0";

   my @ver = (1,16,0);
   my $testcode = <<"---";
#include <mongoc/mongoc.h>
#include <bson/bson.h>

int
main (int argc, char *argv[])
{
#if !MONGOC_CHECK_VERSION($ver[0],$ver[1],$ver[2])
#error libmongoc version check failed
#endif
   mongoc_init ();

   bson_t      child;
   bson_t     *document = bson_new ();
   BSON_APPEND_DOCUMENT_BEGIN (document, "software", &child);
   BSON_APPEND_UTF8 (&child, "name", "polymake");
   bson_append_document_end (document, &child);
   bson_destroy (document);

   mongoc_cleanup ();
   return 0;
}

---

   $add_mongoc_flag = "";
   $add_bson_flag = "";
   $mongoc_tries = 0;
   $bson_tries = 0;
RETRY:
   my $error=Polymake::Configure::build_test_program($testcode, CXXFLAGS => $CXXFLAGS.$add_mongoc_flag.$add_bson_flag, LDFLAGS => $LDFLAGS, LIBS => $LIBS);
   if ($? != 0) {
      if ($error =~ /libmongoc version check failed/) {
         my $verstr = join(".",@ver);
         $options->{polydb} = "" if $usedopts > 0;
         die <<"---";
libmongoc version too old, please install a newer version (>= $verstr) and use
--with-mongoc=prefix to specify the installation directory.
---
      }
      # retry with /usr/include/<subdir> only when no --with flags are used
      if ($usedopts == 0 && $mongoc_tries == 0 && $error =~ /#include.*mongoc\.h.*/ ) {
         $mongoc_tries = 1;
         $add_mongoc_flag .= " -I/usr/include/$lmc1";
         goto RETRY;
      }
      if ($usedopts == 0 && $mongoc_tries == 1 && $error =~ /#include.*mongoc\.h.*/ && ( defined($pkg_base=$Polymake::Configure::BrewBase) or defined($pkg_base=$Polymake::Configure::FinkBase) )) {
         $mongoc_tries = 2;
         $add_mongoc_flag .= " -I$pkg_base/include/$lmc1";
         goto RETRY;
      }
      if ($usedopts == 0 && $bson_tries == 0 && $error =~ /#include.*bson\.h.*/ && $CXXFLAGS !~ /$lb1/) {
         $bson_tries = 1;
         $add_bson_flag .= " -I/usr/include/$lb1";
         goto RETRY;
      }
      if ($usedopts == 0 && $bson_tries == 1 && $error =~ /#include.*bson\.h.*/ && $CXXFLAGS !~ /$lb1/ && ( defined($pkg_base=$Polymake::Configure::BrewBase) or defined($pkg_base=$Polymake::Configure::FinkBase) )) {
         $bson_tries = 2;
         $add_bson_flag .= " -I$pkg_base/include/$lb1";
         goto RETRY;
      }

      $options->{polydb} = "" if $usedopts > 0;
      die <<"---";
Could not compile test program checking for libmongoc and libbson.
Please try using --with-mongoc to specify the installation prefix, or check
./configure --help for more flags to set the include and libary paths.
The complete error log follows:\n\n$error\n
---
   }
   $CXXFLAGS .= $add_mongoc_flag.$add_bson_flag;
   return "CXXFLAGS '$CXXFLAGS', LDFLAGS '$LDFLAGS'";
}
