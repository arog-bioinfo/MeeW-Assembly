"""Derive and validate minimap2 presets for Assembly long-read types."""

import shlex

MINIMAP2_PRESETS = {
    "nano-raw": "map-ont",
    "nano-corr": "map-ont",
    "nano-hq": "map-ont",
    "pacbio-raw": "map-pb",
    "pacbio-corr": "map-pb",
    "pacbio-hifi": "map-hifi",
}

# A preset option can follow other no-argument short options in one token, such
# as ``-axmap-ont``. Options that take an argument consume the rest of a token,
# so a later ``x`` in that argument is not a preset selector.
_NO_ARGUMENT_SHORT_OPTIONS = frozenset("2aSDVXHcQhLYP")


def minimap2_preset(long_read_type):
    """Return the minimap2 mapping preset for an exact long-read type."""
    try:
        return MINIMAP2_PRESETS[long_read_type]
    except (KeyError, TypeError) as error:
        supported = ", ".join(sorted(MINIMAP2_PRESETS))
        raise ValueError(
            f"Unsupported `long_read_type` {long_read_type!r}; expected one of: "
            f"{supported}. Cannot derive a minimap2 preset."
        ) from error


def _shell_arguments(extra):
    try:
        lexer = shlex.shlex(extra, posix=True, punctuation_chars=True)
        lexer.commenters = ""
        lexer.whitespace_split = True
        return list(lexer)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"`minimap2_extra` must be a valid POSIX shell argument string: {error}"
        ) from error


def _preset_option_offset(argument):
    """Return the offset of a minimap2 ``-x`` short option, if present."""
    if len(argument) < 2 or not argument.startswith("-") or argument.startswith("--"):
        return None

    for offset, option in enumerate(argument[1:], start=1):
        if option == "x":
            return offset
        if option not in _NO_ARGUMENT_SHORT_OPTIONS:
            return None
    return None


def validate_minimap2_extra(extra):
    """Reject preset selection while accepting other trusted shell arguments."""
    arguments = _shell_arguments(extra)
    for index, argument in enumerate(arguments):
        offset = _preset_option_offset(argument)
        if offset is None:
            continue

        if argument == "-x" and index == len(arguments) - 1:
            form = "bare `-x` with no value"
        elif argument == "-x":
            form = f"separated preset option `-x {arguments[index + 1]}`"
        elif offset == 1:
            form = f"attached preset option `{argument}`"
        else:
            form = f"bundled preset option `{argument}`"
        raise ValueError(
            f"`minimap2_extra` contains a forbidden {form}. Remove `-x`; the "
            "minimap2 preset is derived automatically from each sample's "
            "`long_read_type`."
        )


def compose_minimap2_extra(long_read_type, extra):
    """Place one derived preset before validated, non-preset extra arguments."""
    preset = minimap2_preset(long_read_type)
    validate_minimap2_extra(extra)
    stripped_extra = extra.strip()
    automatic_preset = f"-x {preset}"
    if not stripped_extra:
        return automatic_preset
    return f"{automatic_preset} {stripped_extra}"
