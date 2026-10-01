//! `mung describe`: the flags front ends build their forms from.

use std::process::Command;

fn mung(args: &[&str], envs: &[(&str, &std::path::Path)]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_mung"))
        .args(args)
        .envs(envs.iter().copied())
        .output()
        .unwrap()
}

fn describe(name: &str) -> std::process::Output {
    mung(&["describe", name], &[])
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
    assert_eq!(arg("gff")["required"], false);
    assert_eq!(arg("species")["value_type"], "text");
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
fn clones_wants_query_unless_from() {
    let out = mung(&["clones", "--out", "x"], &[]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("QUERY"), "{err}");
}

#[test]
fn without_gff_offline_and_uncached_says_what_to_do() {
    let tmp = std::env::temp_dir().join(format!("mung-offline-{}", std::process::id()));
    let (cache, config) = (tmp.join("cache"), tmp.join("config"));
    let out = mung(
        &["infercnv", "--out", "x", "q.zarr.zip"],
        &[
            ("MUNG_OFFLINE", std::path::Path::new("1")),
            ("MUNG_CACHE_DIR", &cache),
            ("MUNG_CONFIG_DIR", &config),
        ],
    );
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("--gff") && err.contains("mung data fetch"),
        "{err}"
    );
}
