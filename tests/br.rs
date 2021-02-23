use std::io::Read;
use std::process::{Command, Stdio};

#[test]
fn count() {
    let mut child = Command::new("./target/debug/br")
        .args(&[
            "-vvv",
            "-i",
            "tests/data/raw.fasta",
            "-o",
            "tests/data/corr.fasta",
            "-k",
            "11",
        ])
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Couldn't create br subprocess");

    if !child.wait().expect("Error durring br run").success() {
        let mut stdout = String::new();
        let mut stderr = String::new();

        child.stdout.unwrap().read_to_string(&mut stdout).unwrap();
        child.stderr.unwrap().read_to_string(&mut stderr).unwrap();

        println!("stdout: {}", stdout);
        println!("stderr: {}", stderr);
        panic!();
    }
}

#[test]
fn solid() {
    let mut child = Command::new("./target/debug/br")
        .args(&[
            "-vvv",
            "-i",
            "tests/data/raw.fasta",
            "-o",
            "tests/data/corr.fasta",
            "-s",
            "tests/data/raw.k11.a2.solid",
        ])
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Couldn't create br subprocess");

    if !child.wait().expect("Error durring br run").success() {
        let mut stdout = String::new();
        let mut stderr = String::new();

        child.stdout.unwrap().read_to_string(&mut stdout).unwrap();
        child.stderr.unwrap().read_to_string(&mut stderr).unwrap();

        println!("stdout: {}", stdout);
        println!("stderr: {}", stderr);
        panic!();
    }
}
