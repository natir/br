/* std use */

/* 3rd party use */

#[cfg(test)]
mod tests {

    #[test]
    fn fasta_count() -> std::io::Result<()> {
        let mut cmd = assert_cmd::Command::cargo_bin("br").unwrap();
        cmd.args(&[
            "-i",
            "tests/data/raw.fasta",
            "-o",
            "tests/data/corr.fasta",
            "fasta",
            "-i",
            "tests/data/raw.fasta",
            "-k",
            "11",
            "first-minimum",
        ]);

        let assert = cmd.assert();

        assert.success().stderr(b"" as &[u8]);

        Ok(())
    }

    #[test]
    fn solid() -> std::io::Result<()> {
        let mut cmd = assert_cmd::Command::cargo_bin("br").unwrap();
        cmd.args(&[
            "-i",
            "tests/data/raw.fasta",
            "-o",
            "tests/data/corr.fasta",
            "solid",
            "-i",
            "tests/data/raw.k11.a2.solid",
            "-f",
            "solid",
        ]);

        let assert = cmd.assert();

        assert.success().stderr(b"" as &[u8]);

        Ok(())
    }

    #[test]
    fn large_kmer() -> std::io::Result<()> {
        let mut cmd = assert_cmd::Command::cargo_bin("br").unwrap();
        cmd.args(&[
            "-i",
            "tests/data/raw.fasta",
            "-o",
            "tests/data/corr.fasta",
            "large-kmer",
            "-i",
            "tests/data/raw.k31.fasta",
            "-f",
            "fasta",
            "-k",
            "31",
        ]);

        let assert = cmd.assert();

        assert.success().stderr(b"" as &[u8]);

        Ok(())
    }
}
