"""seed initial roles and admin

Revision ID: b9b0f7c3d1aa
Revises: 4ce0890619df
Create Date: 2026-04-07 00:00:00.000000

"""
from datetime import datetime
import os

from alembic import op
import sqlalchemy as sa
import bcrypt


# revision identifiers, used by Alembic.
revision = 'b9b0f7c3d1aa'
down_revision = '4ce0890619df'
branch_labels = None
depends_on = None


role_table = sa.table(
    'role',
    sa.column('id', sa.Integer),
    sa.column('name', sa.String),
)

accounts_table = sa.table(
    'accounts',
    sa.column('id', sa.Integer),
    sa.column('username', sa.String),
    sa.column('email', sa.String),
    sa.column('name', sa.String),
    sa.column('surname', sa.String),
    sa.column('password', sa.String),
    sa.column('created_at', sa.DateTime),
    sa.column('updated_at', sa.DateTime),
    sa.column('role_id', sa.Integer),
)


def _admin_seed_values():
    admin_password = os.getenv('NPMINE_WEB_APP_PASSWORD', 'change-me-now')
    admin_email = os.getenv('NPMINE_WEB_APP_EMAIL', 'admin@example.com')
    hashed_password = bcrypt.hashpw(
        admin_password.encode('utf-8'),
        bcrypt.gensalt(),
    ).decode('utf-8')

    return {
        'username': 'admin',
        'email': admin_email,
        'name': 'Admin',
        'surname': 'User',
        'password': hashed_password,
        'created_at': datetime.utcnow(),
        'updated_at': datetime.utcnow(),
        'role_id': 1,
    }


def upgrade():
    bind = op.get_bind()

    existing_role_ids = {
        row[0]
        for row in bind.execute(sa.select(role_table.c.id))
    }

    for role_id, role_name in ((1, 'admin'), (2, 'editor'), (3, 'user')):
        if role_id not in existing_role_ids:
            op.bulk_insert(
                role_table,
                [{'id': role_id, 'name': role_name}],
            )

    admin_exists = bind.execute(
        sa.select(accounts_table.c.id).where(accounts_table.c.username == 'admin')
    ).scalar()

    if not admin_exists:
        op.bulk_insert(accounts_table, [_admin_seed_values()])


def downgrade():
    bind = op.get_bind()

    bind.execute(
        sa.delete(accounts_table).where(accounts_table.c.username == 'admin')
    )
    bind.execute(
        sa.delete(role_table).where(role_table.c.id.in_([1, 2, 3]))
    )
