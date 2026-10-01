//! `mung describe`: the flags front ends build their forms from.

use std::process::Command;

fn describe(name: &str) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_mung"))
        .args(["describe", name])
        .output()
        .unwrap()
}

#[test]
fn describes_the_clones_flags_as_json() {
    let out = describe("clones");
    assert!(out.status.success());
    let v: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["describe"], 1);
    assert_eq!(v["command"], "clones");
    let args = v["args"].as_array().unwrap();
    let arg = |long: &str| args.iter().find(|a| a["long"] == long).unwrap();
    assert_eq!(arg("out")["required"], true);
    assert_eq!(arg("k-max")["default"][0], "8");
    assert_eq!(arg("exclude-chr")["delimiter"], ",");
    assert_eq!(arg("k-max")["value_type"], "unsigned");
    assert_eq!(arg("bin-size")["value_type"], "integer");
    assert_eq!(arg("clip")["value_type"], "number");
    assert_eq!(arg("gff")["value_type"], "text");
    assert_eq!(arg("no-center")["action"], "set_true");
    assert!(arg("engine")["values"]
        .as_array()
        .unwrap()
        .iter()
        .any(|x| x == "bayes"));
    assert!(args
        .iter()
        .any(|a| a["positional"] == true && a["id"] == "query"));
}

#[test]
fn an_unknown_command_is_an_error() {
    assert!(!describe("nope").status.success());
    assert!(!describe("describe").status.success());
}

#[test]
fn clones_wants_gff_and_query_unless_from() {
    let run = |args: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_mung"))
            .arg("clones")
            .args(args)
            .output()
            .unwrap()
    };
    let missing = run(&["--out", "x"]);
    assert!(!missing.status.success());
    let err = String::from_utf8_lossy(&missing.stderr);
    assert!(err.contains("--gff") && err.contains("QUERY"), "{err}");
}
