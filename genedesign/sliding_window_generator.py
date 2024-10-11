def sliding_window_generator(sequence, n_in_scope=3, n_ahead=6, step=None):
    """
    A generator that yields combined in-scope and context-after sequences,
    moving in steps of n_in_scope amino acids.

    :param sequence: The full amino acid sequence (string).
    :param n_in_scope: Number of amino acids in the in-scope segment (default 3).
    :param n_ahead: Number of amino acids after the in-scope segment (default 6).
    :param step: The step size for moving the window (default is n_in_scope).
    :yield: A string representing the combined in-scope and context-after sequence.
    """
    if step is None:
        step = n_in_scope  # Default step size is the size of the in-scope segment

    seq_length = len(sequence)

    for i in range(0, seq_length, step):
        # Extract in-scope amino acids
        in_scope = sequence[i:i + n_in_scope]
        if not in_scope:
            break  # No more amino acids to process

        # Extract context after
        context_after = sequence[i + n_in_scope:i + n_in_scope + n_ahead]

        # Combine in-scope and context after into a single string
        window = in_scope + context_after

        # Yield the window string
        yield window