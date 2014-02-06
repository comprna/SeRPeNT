# unit_test.pl
# Parses the cases file and executes the cases reporting FAILURE or SUCCESS
#
# Parameters:
#   $ARGV[0] : Path of the cases file
#   $ARGV[1] : Path of the report file

my ($cases_file, $report_file, @arguments) = @ARGV;
my $result = 1;

open(REPORT, " > $report_file");

open(CASES, " < $cases_file");
while(<CASES>) {
  chomp;
  next if (m/^#/);
  my ($key, $value) = split "\t";

  print REPORT "$value\n" if($key eq "UTM");
  print REPORT "\t$value\t" if ($key eq "CAE");
  if ($key eq "CMD") {
    for (my $i = 1; $i <= scalar(@arguments); $i++) {
      $value =~ s/\$$i/$arguments[$i - 1]/g;
    }
    $result = `$value`;
    my @lines = split("\n", $result);
    $result = $lines[-1];
  }
  if ($key eq "OUT") {
    print REPORT "[FAIL]\n" if ($result ne $value);
    print REPORT "[PASS]\n" if ($result eq $value);
  }
}
close(CASES);
close(REPORT);
