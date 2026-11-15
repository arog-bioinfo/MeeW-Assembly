"""Shared derivation and validation for Assembly genome context."""

import shlex

ALLOWED_GENOME_CONTEXTS = ("metagenome", "isolate")
DEFAULT_IMPORTED_GENOME_CONTEXT = "metagenome"


def resolve_genome_context(value):
    """Validate a supplied genome context."""
    if value not in ALLOWED_GENOME_CONTEXTS:
        allowed = " or ".join(f"`{context}`" for context in ALLOWED_GENOME_CONTEXTS)
        raise ValueError(f"`genome_context` must be exactly {allowed}; got {value!r}.")
    return value


def resolve_config_genome_context(config):
    """Fallback only when an imported config omits genome_context entirely."""
    if "genome_context" not in config:
        return DEFAULT_IMPORTED_GENOME_CONTEXT
    return resolve_genome_context(config["genome_context"])


def quast_tool(genome_context):
    """Select QUAST's executable exclusively from genome context."""
    return (
        "metaquast"
        if resolve_genome_context(genome_context) == "metagenome"
        else "quast"
    )


def reject_meta_extra(extra, config_key):
    """Reject a duplicate/context-conflicting --meta in a free-form extra string."""
    try:
        lexer = shlex.shlex(extra, posix=True, punctuation_chars=True)
        lexer.commenters = ""
        lexer.whitespace_split = True
        arguments = list(lexer)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"`{config_key}` must be a valid shell argument string: {error}"
        ) from error

    if any(
        argument == "--meta" or argument.startswith("--meta=") for argument in arguments
    ):
        raise ValueError(
            f"Do not set `--meta` in `{config_key}`; set `genome_context: metagenome` "
            "to enable metagenome assembly mode."
        )
