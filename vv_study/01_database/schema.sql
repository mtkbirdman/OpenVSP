-- VSPAERO 3.50.5 V&V database specification DDL v1.0
-- SQLite 3, enable foreign keys on EVERY connection. Apply schema in a transaction.
PRAGMA foreign_keys = ON;

CREATE TABLE study (
    study_id TEXT PRIMARY KEY,
    name TEXT NOT NULL,
    plan_version TEXT NOT NULL,
    plan_sha256 TEXT NOT NULL CHECK(length(plan_sha256)=64),
    target_openvsp_version TEXT NOT NULL CHECK(target_openvsp_version='3.50.5'),
    config_sha256 TEXT NOT NULL CHECK(length(config_sha256)=64),
    schema_version TEXT NOT NULL,
    created_utc TEXT NOT NULL
);
CREATE TABLE test_definition (
    test_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    title TEXT NOT NULL,
    requirement TEXT NOT NULL CHECK(requirement IN ('REQUIRED','CONDITIONAL')),
    execution_mode TEXT NOT NULL CHECK(execution_mode IN ('AUTO','SEMI_AUTO','MANUAL')),
    category TEXT NOT NULL,
    description TEXT NOT NULL
);
CREATE TABLE geometry (
    geometry_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    stage TEXT NOT NULL,
    physical_model_json TEXT NOT NULL CHECK(json_valid(physical_model_json)),
    representation_json TEXT NOT NULL CHECK(json_valid(representation_json)),
    source_vsp3_sha256 TEXT NOT NULL CHECK(length(source_vsp3_sha256)=64),
    geometry_equivalence_status TEXT NOT NULL DEFAULT 'UNVERIFIED'
      CHECK(geometry_equivalence_status IN ('UNVERIFIED','REVIEWED','NOT_EQUIVALENT')),
    remarks TEXT
);
CREATE TABLE flight_condition (
    condition_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    mach REAL NOT NULL CHECK(mach>=0),
    re_cref REAL NOT NULL CHECK(re_cref>0),
    alpha_deg REAL NOT NULL,
    beta_deg REAL NOT NULL,
    controls_deg_json TEXT NOT NULL CHECK(json_valid(controls_deg_json)),
    trim_mode TEXT NOT NULL CHECK(trim_mode IN ('FIXED_CONTROL','TRIMMED')),
    sref REAL NOT NULL CHECK(sref>0),
    bref REAL NOT NULL CHECK(bref>0),
    cref REAL NOT NULL CHECK(cref>0),
    moment_x REAL NOT NULL,
    moment_y REAL NOT NULL,
    moment_z REAL NOT NULL,
    reference_length_unit TEXT NOT NULL,
    axes_convention_id TEXT NOT NULL,
    condition_hash TEXT NOT NULL CHECK(length(condition_hash)=64)
);
CREATE TABLE numerical_setting (
    setting_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    analysis_mode TEXT NOT NULL,
    mesh_level TEXT NOT NULL,
    requested_json TEXT NOT NULL CHECK(json_valid(requested_json)),
    setting_hash TEXT NOT NULL CHECK(length(setting_hash)=64),
    notes TEXT
);
CREATE TABLE planned_case (
    case_id TEXT PRIMARY KEY,
    test_id TEXT NOT NULL REFERENCES test_definition(test_id),
    role TEXT NOT NULL,
    geometry_id TEXT REFERENCES geometry(geometry_id),
    condition_id TEXT REFERENCES flight_condition(condition_id),
    setting_id TEXT REFERENCES numerical_setting(setting_id),
    execute_policy TEXT NOT NULL CHECK(execute_policy IN ('EXECUTE','REUSE','IMPORT','MANUAL','UNRESOLVED')),
    expected_outputs_json TEXT NOT NULL CHECK(json_valid(expected_outputs_json)),
    UNIQUE(test_id,role)
);
CREATE TABLE run (
    run_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    geometry_id TEXT NOT NULL REFERENCES geometry(geometry_id),
    condition_id TEXT NOT NULL REFERENCES flight_condition(condition_id),
    setting_id TEXT NOT NULL REFERENCES numerical_setting(setting_id),
    status TEXT NOT NULL CHECK(status IN (
      'PLANNED','PRECHECKED','RUNNING','CAPTURED','PARSED','SUCCEEDED',
      'PRECHECK_FAILED','EXECUTION_FAILED','OUTPUT_INVALID','INTERRUPTED','CANCELLED')),
    run_dir TEXT NOT NULL UNIQUE,
    fingerprint_sha256 TEXT NOT NULL CHECK(length(fingerprint_sha256)=64),
    binary_sha256 TEXT NOT NULL CHECK(length(binary_sha256)=64),
    runner_sha256 TEXT NOT NULL CHECK(length(runner_sha256)=64),
    requested_snapshot_json TEXT NOT NULL CHECK(json_valid(requested_snapshot_json)),
    effective_snapshot_json TEXT CHECK(effective_snapshot_json IS NULL OR json_valid(effective_snapshot_json)),
    generated_snapshot_json TEXT CHECK(generated_snapshot_json IS NULL OR json_valid(generated_snapshot_json)),
    environment_json TEXT NOT NULL CHECK(json_valid(environment_json)),
    exit_code INTEGER,
    error_code TEXT,
    created_utc TEXT NOT NULL,
    started_utc TEXT,
    ended_utc TEXT,
    retry_of_run_id TEXT REFERENCES run(run_id)
);
CREATE TABLE case_point (
    point_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES run(run_id),
    ordinal INTEGER NOT NULL CHECK(ordinal>=0),
    target_condition_id TEXT NOT NULL REFERENCES flight_condition(condition_id),
    observed_condition_json TEXT NOT NULL CHECK(json_valid(observed_condition_json)),
    point_status TEXT NOT NULL CHECK(point_status IN ('MATCHED','MISMATCH','MISSING','UNVERIFIED')),
    UNIQUE(run_id,ordinal)
);
CREATE TABLE case_run_link (
    case_id TEXT NOT NULL REFERENCES planned_case(case_id),
    point_id TEXT NOT NULL REFERENCES case_point(point_id),
    linked_utc TEXT NOT NULL,
    link_reason TEXT NOT NULL,
    PRIMARY KEY(case_id,point_id)
);
CREATE TABLE dataset (
    dataset_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    kind TEXT NOT NULL CHECK(kind IN
      ('VSPAERO','PANEL','NS','WIND_TUNNEL','THEORY','EXTERNAL_CD0','DERIVED')),
    point_id TEXT REFERENCES case_point(point_id),
    condition_id TEXT REFERENCES flight_condition(condition_id),
    geometry_id TEXT REFERENCES geometry(geometry_id),
    usage_class TEXT NOT NULL DEFAULT 'REFERENCE'
      CHECK(usage_class IN ('CALIBRATION','VALIDATION','REFERENCE','NOT_APPLICABLE')),
    evidence_class TEXT CHECK(evidence_class IS NULL OR evidence_class IN ('A','B','C')),
    representation TEXT NOT NULL,
    source_description TEXT NOT NULL,
    metadata_json TEXT NOT NULL CHECK(json_valid(metadata_json)),
    method_id TEXT NOT NULL,
    method_version TEXT NOT NULL,
    created_utc TEXT NOT NULL,
    CHECK ((kind='VSPAERO' AND point_id IS NOT NULL)
        OR (kind!='VSPAERO'))
);
CREATE TABLE artifact (
    artifact_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    run_id TEXT REFERENCES run(run_id),
    dataset_id TEXT REFERENCES dataset(dataset_id),
    role TEXT NOT NULL,
    relative_path TEXT NOT NULL UNIQUE,
    sha256 TEXT NOT NULL CHECK(length(sha256)=64),
    size_bytes INTEGER NOT NULL CHECK(size_bytes>=0),
    media_type TEXT,
    created_utc TEXT NOT NULL,
    integrity_status TEXT NOT NULL DEFAULT 'VERIFIED'
      CHECK(integrity_status IN ('VERIFIED','MISSING','CHANGED'))
);
CREATE TABLE quantity (
    quantity_id TEXT PRIMARY KEY,
    dataset_id TEXT NOT NULL REFERENCES dataset(dataset_id),
    qoi TEXT NOT NULL,
    component TEXT NOT NULL DEFAULT 'TOTAL',
    variant TEXT NOT NULL DEFAULT 'DEFAULT',
    native_name TEXT NOT NULL,
    native_value REAL,
    native_unit TEXT,
    canonical_value REAL,
    canonical_unit TEXT,
    convention_id TEXT,
    normalization_method_id TEXT,
    source_artifact_id TEXT REFERENCES artifact(artifact_id),
    source_locator TEXT,
    value_status TEXT NOT NULL CHECK(value_status IN
      ('AVAILABLE','NATIVE_ONLY','NOT_CALCULATED','NOT_AVAILABLE','NOT_APPLICABLE','INVALID')),
    CHECK ((value_status='AVAILABLE' AND canonical_value IS NOT NULL)
        OR (value_status!='AVAILABLE')),
    UNIQUE(dataset_id,qoi,component,variant)
);
CREATE TABLE quantity_dependency (
    derived_quantity_id TEXT NOT NULL REFERENCES quantity(quantity_id),
    input_quantity_id TEXT NOT NULL REFERENCES quantity(quantity_id),
    role TEXT NOT NULL,
    PRIMARY KEY(derived_quantity_id,input_quantity_id,role),
    CHECK(derived_quantity_id<>input_quantity_id)
);
CREATE TABLE acceptance_rule (
    rule_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    qoi TEXT NOT NULL,
    level TEXT NOT NULL CHECK(level IN ('NUMERICAL','ENGINEERING','APPLICABILITY')),
    policy_json TEXT NOT NULL CHECK(json_valid(policy_json)),
    approval TEXT NOT NULL CHECK(approval IN ('PROVISIONAL','APPROVED','RETIRED')),
    rationale TEXT NOT NULL,
    reviewer TEXT,
    approved_utc TEXT
);
CREATE TABLE comparison (
    comparison_id TEXT PRIMARY KEY,
    test_id TEXT NOT NULL REFERENCES test_definition(test_id),
    predicted_quantity_id TEXT NOT NULL REFERENCES quantity(quantity_id),
    reference_quantity_id TEXT NOT NULL REFERENCES quantity(quantity_id),
    comparison_kind TEXT NOT NULL CHECK(comparison_kind IN
      ('DIRECT','CONVERTED','INTERPOLATED','ILLUSTRATIVE')),
    delta_value REAL,
    relative_delta REAL,
    standard_uncertainty REAL,
    expanded_uncertainty REAL,
    uncertainty_metadata_json TEXT NOT NULL CHECK(json_valid(uncertainty_metadata_json)),
    rule_id TEXT REFERENCES acceptance_rule(rule_id),
    verdict TEXT NOT NULL CHECK(verdict IN
      ('PASS','FAIL','INCONCLUSIVE','NOT_EVALUATED')),
    rationale TEXT NOT NULL,
    created_utc TEXT NOT NULL
);
CREATE TABLE assessment (
    assessment_id TEXT PRIMARY KEY,
    study_id TEXT NOT NULL REFERENCES study(study_id),
    test_id TEXT REFERENCES test_definition(test_id),
    case_id TEXT REFERENCES planned_case(case_id),
    qoi TEXT NOT NULL,
    applicability_domain_json TEXT NOT NULL CHECK(json_valid(applicability_domain_json)),
    numerical_status TEXT NOT NULL CHECK(numerical_status IN
      ('PASS','FAIL','INCONCLUSIVE','NOT_EVALUATED')),
    engineering_status TEXT NOT NULL CHECK(engineering_status IN
      ('PASS','FAIL','INCONCLUSIVE','NOT_EVALUATED')),
    disposition TEXT NOT NULL CHECK(disposition IN
      ('APPLICABLE','CONDITIONAL','NOT_APPLICABLE','INSUFFICIENT_EVIDENCE')),
    rule_ids_json TEXT NOT NULL CHECK(json_valid(rule_ids_json)),
    evidence_ids_json TEXT NOT NULL CHECK(json_valid(evidence_ids_json)),
    restriction_json TEXT NOT NULL CHECK(json_valid(restriction_json)),
    rationale TEXT NOT NULL,
    reviewer TEXT,
    approved_utc TEXT,
    created_utc TEXT NOT NULL
);
CREATE TABLE run_event (
    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id TEXT NOT NULL REFERENCES run(run_id),
    occurred_utc TEXT NOT NULL,
    previous_status TEXT,
    next_status TEXT NOT NULL,
    event_type TEXT NOT NULL,
    event_code TEXT NOT NULL,
    message TEXT NOT NULL,
    details_json TEXT NOT NULL CHECK(json_valid(details_json))
);
CREATE INDEX idx_cases_test ON planned_case(test_id);
CREATE INDEX idx_runs_input ON run(geometry_id,condition_id,setting_id,status);
CREATE INDEX idx_points_run ON case_point(run_id);
CREATE INDEX idx_case_run_point ON case_run_link(point_id);
CREATE INDEX idx_dataset_point ON dataset(point_id,kind);
CREATE INDEX idx_quantity_qoi ON quantity(qoi,dataset_id);
CREATE INDEX idx_comparison_test ON comparison(test_id,verdict);
CREATE INDEX idx_assessment_scope ON assessment(test_id,qoi,disposition);
CREATE INDEX idx_artifact_run ON artifact(run_id);
CREATE INDEX idx_run_event ON run_event(run_id,occurred_utc);

-- Additive migration v1: frozen evaluation snapshots and schema migrations.
CREATE TABLE IF NOT EXISTS evaluation_selection (
 label TEXT PRIMARY KEY,
 test_id TEXT NOT NULL REFERENCES test_definition(test_id),
 selection_json TEXT NOT NULL CHECK(json_valid(selection_json)),
 selection_sha256 TEXT NOT NULL CHECK(length(selection_sha256)=64),
 created_utc TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS schema_migration (
 migration_id TEXT PRIMARY KEY,
 applied_utc TEXT NOT NULL,
 notes TEXT NOT NULL
);
INSERT OR IGNORE INTO schema_migration VALUES('EXTENSION_1', strftime('%Y-%m-%dT%H:%M:%fZ','now'), 'Immutable evaluation selection; no raw records changed');
