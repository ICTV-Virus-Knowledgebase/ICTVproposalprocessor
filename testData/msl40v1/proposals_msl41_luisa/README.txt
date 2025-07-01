This file has a spurious value in cell AA1048536, proposed_taxonomy.subfamily, but no value in the "Change" column (or any other for that matter.)

That caused a crash. Now fixed to emit an error, and skip processing that row.
