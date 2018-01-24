package FSSP;

# global variables to be visible by outsiders

use strict;
use DBI;

my $dbname= 'dali';
my $default_socket='/var/lib/mysql/mysql.sock'; # /tmp/mysql-5.0.sock
my $default_host  = 'localhost';
my $default_port  = 3306;

my $dbh;
my $sth;                # contains the last statement handle
my $rv;                 # contains the last result value

##################################################
## Subroutines used for connecting to the database
##################################################
sub Connect {
  my ($user, $password, $host, $socket, $port) = @_;

  if (!$user) {
    $user       = 'wwwrun';
    $password   = 'wwwrun';
  }

  # $ENV{'MYSQL_UNIX_PORT'} = $socket; ???

  #  Connect to DBI
  my $dbs = "DBI:mysql:$dbname";
  if ($host) { $dbs .= ";host=$host"; } else { if ($default_host) { $dbs .= ";host=$default_host";}}
  if ($socket) { $dbs .= ";mysql_socket=$socket"; } else { if ($default_socket) { $dbs .= ";mysql_socket=$default_socket";}}
  if ($port) { $dbs .= ";port=$port";} else { if ($default_port) { $dbs .= ";port=$default_port";}}

  return ($dbh = DBI->connect($dbs, $user, $password));
}

sub DisConnect {

#  Disconnect from DBI
  return $dbh->disconnect();
}

##############################################################

####################################################
# Title:        Execute()
# Function:     Execute an SQL command
####################################################
sub Execute {
  my $command = shift;

  $sth = $dbh->prepare("$command");
  return ($rv = $sth->execute());
}

####################################################
# Title:        Select()
# Function:     Execute an SQL command and return a result structure
####################################################
sub Select {
  my $command = shift;

  $sth = $dbh->prepare("$command");
  $rv = $sth->execute();
  if ($rv) {
    return BuildResultStructure();
  } else {
    return 0;
  }
}

sub CheckForSpecialCharacters {
  my $string = shift;
  $string =~ s/([\'\"])/\\\\\\$1/;
  return $string;
}

sub Quote {
  my $string = shift;
  return $dbh->quote($string);
}

####################################################################
# remarks:
# fetchrow_hashref is used for better readability
# although it is not supposed to be as fast as
# fetchrow_arrayref. Since this part of Picasso is
# used for display, there shouldn't be that many rows
# anyway. For high-performance data processing other
# routines should be written.

sub BuildResultStructure {
  my @result;
  my ($hash_ref);

  while($hash_ref = $sth->fetchrow_hashref) {
    push @result, $hash_ref;
  }

  return \@result; # reference to an array of hashreferences
}

#########################################################################
sub JoinResults {
  my $firstquery = shift;
  my $query;
  my $entry;

  foreach $query (@_) {
    foreach $entry (@$query) {
      push @$firstquery, $entry;
    }
  }

  return $firstquery;
}

##########################################################
### Subroutines used for specific queries: used by tools
##########################################################

sub Get_Sequence {
  my ( $cd1 ) = @_;
  Execute("SELECT sequence FROM sequence WHERE cd1 = \'$cd1\'");
  return BuildResultStructure();
}

sub Get_Alignment {
  my ( $cd1,$cd2 ) = @_;
  Execute("SELECT cd1_blocks,cd2_blocks FROM pairs WHERE cd1 = \'$cd1\' AND cd2 = \'$cd2\'");
  return BuildResultStructure();
}

return 1;
