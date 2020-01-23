//! Collection of (V)OPRF specific errors

use std::io::{Error, ErrorKind};

/// Unsupported curve type specified
pub fn unsupported() -> Error { Error::new(ErrorKind::Other, "Curve choice is not supported") }

/// Error deserializing bytes into a valid group element object
pub fn deserialization() -> Error { Error::new(ErrorKind::Other, "Failed to deserialize") }

/// Indicates that an internal error occurred
pub fn internal() -> Error { Error::new(ErrorKind::Other, "Internal error occurred") }
