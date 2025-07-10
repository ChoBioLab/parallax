// modules/local/channel_monitor.nf
process CHANNEL_MONITOR {
    tag "monitor"
    label 'process_single'

    input:
    val(stage_name)
    val(channel_count)
    val(expected_count)

    output:
    path "channel_report_${stage_name}.json", emit: report

    script:
    """
    #!/usr/bin/env python3
    import json
    from datetime import datetime

    report = {
        'stage': '${stage_name}',
        'timestamp': datetime.now().isoformat(),
        'channel_count': ${channel_count},
        'expected_count': ${expected_count},
        'status': 'pass' if ${channel_count} == ${expected_count} else 'fail'
    }

    if report['status'] == 'fail':
        print(f"WARNING: Channel count mismatch at {report['stage']}")
        print(f"Expected: {report['expected_count']}, Got: {report['channel_count']}")

    with open('channel_report_${stage_name}.json', 'w') as f:
        json.dump(report, f, indent=2)
    """
}

